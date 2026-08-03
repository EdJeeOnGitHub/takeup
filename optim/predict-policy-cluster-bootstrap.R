#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")

parameter_csv <- policy_option_value(args, "--parameter-csv")
distance_data <- policy_option_value(args, "--distance-data", "optim/data/full-many-pots-experiment.rds")
output_path <- policy_option_value(args, "--output-path")
distance_cap <- as.numeric(policy_option_value(args, "--distance-cap", "3500"))
num_cores <- as.integer(policy_option_value(args, "--num-cores", "1"))
num_replicates <- as.integer(policy_option_value(args, "--num-replicates", "210"))
if (is.null(parameter_csv) || is.null(output_path)) {
  stop("--parameter-csv and --output-path are required.", call. = FALSE)
}

parameters <- read.csv(parameter_csv, stringsAsFactors = FALSE)
parameters <- parameters[seq_len(min(num_replicates, nrow(parameters))), ]
required <- c(
  "draw", "replicate", "beta_control", "dist_beta", "mu_control",
  "mu_bracelet", "mu_dist_control", "mu_dist_bracelet", "base_mu_rep",
  "u_sd", "total_error_sd"
)
if (!all(required %in% names(parameters)) || any(!is.finite(as.matrix(parameters[required])))) {
  stop("Policy parameter input is incomplete or non-finite.", call. = FALSE)
}

distance_object <- readRDS(distance_data)
edges <- distance_object$long_distance_mat
edges <- edges[edges$dist <= distance_cap, c("index_i", "index_j", "dist", "dist_km")]
names(edges) <- c("village_i", "pot_j", "distance", "distance_km")
edges <- edges[order(edges$village_i, edges$pot_j), ]
if (!setequal(unique(edges$village_i), distance_object$village_df$id)) {
  stop("At least one village has no feasible PoT.", call. = FALSE)
}
unique_distances <- sort(unique(edges$distance))
parameters$sd_of_dist <- distance_object$sd_of_dist

started <- Sys.time()
draw_predictions <- parallel::mclapply(
  seq_len(nrow(parameters)),
  function(index) predict_policy_draw(parameters[index, ], unique_distances),
  mc.cores = num_cores,
  mc.preschedule = TRUE
)
curves <- do.call(rbind, draw_predictions)
if (any(!is.finite(curves$demand)) || any(curves$demand < 0 | curves$demand > 1)) {
  stop("Non-finite or invalid predicted demand.", call. = FALSE)
}

# Evaluate the experimental Control allocation at its exact assigned distances.
experimental_distances <- distance_object$village_df$dist.to.pot
experimental <- do.call(rbind, parallel::mclapply(
  seq_len(nrow(parameters)),
  function(index) {
    prediction <- predict_policy_draw(
      parameters[index, ], experimental_distances,
      policy_scenarios[policy_scenarios$scenario == "control", ]
    )
    prediction$village_i <- distance_object$village_df$id
    prediction
  },
  mc.cores = num_cores,
  mc.preschedule = TRUE
))

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
saveRDS(edges, file.path(output_path, "policy-feasible-edges.rds"), compress = FALSE)
saveRDS(curves, file.path(output_path, "policy-demand-curves.rds"), compress = FALSE)
saveRDS(experimental, file.path(output_path, "policy-experimental-demand.rds"), compress = FALSE)
write.csv(data.frame(
  draws = nrow(parameters),
  scenarios = nrow(policy_scenarios),
  feasible_edges = nrow(edges),
  unique_distances = length(unique_distances),
  dense_rows_avoided = nrow(distance_object$long_distance_mat) * nrow(parameters) * nrow(policy_scenarios) - nrow(curves),
  fixedpoint_fallbacks = sum(vapply(
    split(curves$fixedpoint_fallbacks, interaction(curves$draw, curves$scenario_id)),
    function(value) unique(value)[1L], numeric(1)
  )),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
), file.path(output_path, "policy-prediction-status.csv"), row.names = FALSE)
message("Wrote exact sparse policy demand for ", nrow(parameters), " draws.")
