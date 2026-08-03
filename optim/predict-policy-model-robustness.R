#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")

parameter_rds <- policy_option_value(args, "--parameter-rds")
distance_data <- policy_option_value(args, "--distance-data", "optim/data/full-many-pots-experiment.rds")
output_path <- policy_option_value(args, "--output-path")
distance_cap <- as.numeric(policy_option_value(args, "--distance-cap", "3500"))
num_cores <- as.integer(policy_option_value(args, "--num-cores", "1"))
max_draws <- as.integer(policy_option_value(args, "--max-draws", "0"))
if (is.null(parameter_rds) || is.null(output_path)) {
  stop("--parameter-rds and --output-path are required.", call. = FALSE)
}

parameter_object <- readRDS(parameter_rds)
parameters <- parameter_object$parameters
if (max_draws > 0L) parameters <- parameters[seq_len(min(max_draws, nrow(parameters))), ]
cluster_shock <- parameter_object$cluster_shock
if (!is.null(cluster_shock)) cluster_shock <- cluster_shock[seq_len(nrow(parameters)), , drop = FALSE]

distance_object <- readRDS(distance_data)
if (!is.null(cluster_shock)) {
  cluster_index <- match(
    distance_object$village_df$cluster.id,
    parameter_object$cluster_external_id
  )
  if (anyNA(cluster_index) || anyDuplicated(cluster_index) ||
      length(cluster_index) != ncol(cluster_shock)) {
    stop("Policy villages do not align with the Stan cluster index map.", call. = FALSE)
  }
  cluster_shock <- cluster_shock[, cluster_index, drop = FALSE]
  cluster_crosswalk <- data.frame(
    policy_village_i = distance_object$village_df$id,
    cluster.id = distance_object$village_df$cluster.id,
    stan_cluster_id = cluster_index
  )
} else {
  cluster_crosswalk <- NULL
}
edges <- distance_object$long_distance_mat
edges <- edges[edges$dist <= distance_cap, c("index_i", "index_j", "dist", "dist_km")]
names(edges) <- c("village_i", "pot_j", "distance", "distance_km")
edges <- edges[order(edges$village_i, edges$pot_j), ]
parameters$sd_of_dist <- distance_object$sd_of_dist

started <- Sys.time()
draw_predictions <- parallel::mclapply(seq_len(nrow(parameters)), function(index) {
  tryCatch({
    parameter <- as.list(parameters[index, ])
    prediction <- if (!is.null(cluster_shock)) {
      parameter$cluster_shock <- cluster_shock[index, ]
      predict_policy_draw(parameter, edges$distance, village_ids = edges$village_i)
    } else {
      predict_policy_draw(parameter, sort(unique(edges$distance)))
    }
    if (is.null(cluster_shock)) prediction else prediction$demand
  }, error = function(error) structure(
    list(index = index, message = conditionMessage(error)),
    class = "policy_prediction_error"
  ))
}, mc.cores = num_cores, mc.preschedule = TRUE)
prediction_errors <- vapply(
  draw_predictions, inherits, logical(1), what = "policy_prediction_error"
)
if (any(prediction_errors)) {
  details <- vapply(draw_predictions[prediction_errors], function(error) {
    paste0("draw ", error$index, ": ", error$message)
  }, character(1))
  stop("Policy prediction failed: ", paste(details, collapse = "; "), call. = FALSE)
}
if (is.null(cluster_shock)) {
  curves <- do.call(rbind, draw_predictions)
  if (any(!is.finite(curves$demand)) || any(curves$demand < 0 | curves$demand > 1)) {
    stop("Non-finite or invalid predicted demand.", call. = FALSE)
  }
  edge_demand <- NULL
} else {
  edge_demand <- do.call(rbind, draw_predictions)
  expected_columns <- nrow(edges) * nrow(policy_scenarios)
  finite_edge_demand <- edge_demand[is.finite(edge_demand)]
  if (ncol(edge_demand) != expected_columns ||
      any(finite_edge_demand < 0 | finite_edge_demand > 1)) {
    stop("Invalid compact edge-demand matrix.", call. = FALSE)
  }
  curves <- NULL
}

experimental_distances <- distance_object$village_df$dist.to.pot
experimental <- do.call(rbind, parallel::mclapply(seq_len(nrow(parameters)), function(index) {
  parameter <- as.list(parameters[index, ])
  if (!is.null(cluster_shock)) parameter$cluster_shock <- cluster_shock[index, ]
  prediction <- predict_policy_draw(
    parameter, experimental_distances,
    policy_scenarios[policy_scenarios$scenario == "control", ],
    village_ids = if (is.null(cluster_shock)) NULL else distance_object$village_df$id
  )
  prediction$village_i <- distance_object$village_df$id
  prediction
}, mc.cores = num_cores, mc.preschedule = TRUE))

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
if (!is.null(cluster_crosswalk)) {
  write.csv(
    cluster_crosswalk, file.path(output_path, "policy-cluster-crosswalk.csv"),
    row.names = FALSE
  )
}
saveRDS(edges, file.path(output_path, "policy-feasible-edges.rds"), compress = FALSE)
if (is.null(edge_demand)) {
  saveRDS(curves, file.path(output_path, "policy-demand-curves.rds"), compress = FALSE)
} else {
  saveRDS(
    edge_demand, file.path(output_path, "policy-edge-demand-matrix.rds"),
    compress = FALSE
  )
  write.csv(
    parameters[, c("draw", "replicate", "chain", "iteration", "model_id"), drop = FALSE],
    file.path(output_path, "policy-edge-demand-draw-map.csv"), row.names = FALSE
  )
}
saveRDS(experimental, file.path(output_path, "policy-experimental-demand.rds"), compress = FALSE)
write.csv(data.frame(
  model_id = parameters$model_id[1], draws = nrow(parameters),
  scenarios = nrow(policy_scenarios), feasible_edges = nrow(edges),
  edge_specific = !is.null(cluster_shock),
  prediction_values = if (is.null(edge_demand)) nrow(curves) else length(edge_demand),
  undefined_prediction_values = if (is.null(edge_demand)) 0L else sum(!is.finite(edge_demand)),
  storage_mb = if (is.null(edge_demand)) as.numeric(object.size(curves)) / 1024^2 else
    as.numeric(object.size(edge_demand)) / 1024^2,
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
), file.path(output_path, "policy-prediction-status.csv"), row.names = FALSE)
message("Wrote policy demand for ", nrow(parameters), " draws of ", parameters$model_label[1], ".")
