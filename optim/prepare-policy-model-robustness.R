#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")

input_csv <- policy_option_value(args, "--input-csv")
output_path <- policy_option_value(args, "--output-path")
model_id <- policy_option_value(args, "--model-id")
model_label <- policy_option_value(args, "--model-label", model_id)
model_family <- policy_option_value(args, "--model-family", "gaussian")
lambda_structure <- policy_option_value(args, "--lambda-structure", "common")
lambda_prior_sd <- as.numeric(policy_option_value(args, "--lambda-prior-sd", "0.25"))
workspace <- policy_option_value(
  args, "--workspace", "data/stan_analysis_data/dist_fit106.RData"
)
if (any(vapply(list(input_csv, output_path, model_id), is.null, logical(1)))) {
  stop("--input-csv, --output-path, and --model-id are required.", call. = FALSE)
}

draws <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(draws) == 0L) stop("Compact draw CSV contains no draws.", call. = FALSE)
required_metadata <- c("chain", "iteration", "source_csv")
if (!all(required_metadata %in% names(draws))) {
  stop("Compact draw CSV lacks chain/iteration/source_csv.", call. = FALSE)
}

parameters <- do.call(rbind, lapply(seq_len(nrow(draws)), function(index) {
  canonical_policy_draw(
    draws[index, , drop = FALSE], draw_id = index,
    chain = draws$chain[index], model_id = model_id, model_label = model_label,
    model_family = model_family, lambda_structure = lambda_structure,
    lambda_log_ratio_sd_prior = lambda_prior_sd,
    source_csv = draws$source_csv[index]
  )
}))
parameters$iteration <- draws$iteration

shock_columns <- grep("^core_cluster_shock_raw[.]", names(draws), value = TRUE)
cluster_shock <- NULL
cluster_external_id <- NULL
if (length(shock_columns)) {
  shock_columns <- shock_columns[order(as.integer(sub(".*[.]", "", shock_columns)))]
  shock_sd <- as.numeric(draws[["core_cluster_shock_sd.1"]])
  cluster_shock <- as.matrix(draws[, shock_columns, drop = FALSE]) * shock_sd
  if (any(!is.finite(cluster_shock))) stop("Invalid cluster shock draws.", call. = FALSE)
  if (!file.exists(workspace)) stop("Cluster mapping workspace not found.", call. = FALSE)
  fit_environment <- new.env(parent = emptyenv())
  suppressWarnings(load(workspace, envir = fit_environment))
  analysis_data <- fit_environment$stan_data$analysis_data
  cluster_map <- unique(analysis_data[, c("cluster_id", "cluster.id")])
  cluster_map <- cluster_map[order(cluster_map$cluster_id), ]
  if (!identical(cluster_map$cluster_id, seq_len(ncol(cluster_shock))) ||
      anyDuplicated(cluster_map$cluster.id)) {
    stop("Stan cluster index mapping is incomplete or ambiguous.", call. = FALSE)
  }
  cluster_external_id <- cluster_map$cluster.id
}

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    parameters = parameters, cluster_shock = cluster_shock,
    cluster_external_id = cluster_external_id
  ),
  file.path(output_path, "policy-model-parameters.rds"), compress = FALSE
)
write.csv(parameters, file.path(output_path, "policy-model-parameters.csv"), row.names = FALSE)
if (!is.null(cluster_external_id)) {
  write.csv(
    data.frame(
      stan_cluster_id = seq_along(cluster_external_id),
      cluster.id = cluster_external_id
    ),
    file.path(output_path, "policy-stan-cluster-map.csv"), row.names = FALSE
  )
}
write.csv(data.frame(
  model_id = model_id, model_label = model_label, model_family = model_family,
  lambda_structure = lambda_structure, lambda_prior_sd = lambda_prior_sd,
  draws = nrow(parameters), chains = length(unique(parameters$chain)),
  cluster_shocks = if (is.null(cluster_shock)) 0L else ncol(cluster_shock)
), file.path(output_path, "policy-model-parameter-status.csv"), row.names = FALSE)
message("Prepared ", nrow(parameters), " policy draws for ", model_label, ".")
