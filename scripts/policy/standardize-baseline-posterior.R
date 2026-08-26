#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")
input_path <- policy_option_value(args, "--input-path", "build/policy/posterior")
output_path <- policy_option_value(args, "--output-path")
if (is.null(output_path)) stop("--output-path is required.", call. = FALSE)

required <- c(
  "policy-posterior-parameters.csv", "policy-edge-demand-draw-map.csv",
  "policy-edge-demand-matrix.rds", "policy-experimental-demand.rds",
  "policy-feasible-edges.rds", "policy-prediction-status.csv"
)
missing <- required[!file.exists(file.path(input_path, required))]
if (length(missing)) stop("Missing baseline files: ", paste(missing, collapse = ", "))
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(input_path, required[-1L]), output_path, overwrite = TRUE)
parameters <- read.csv(file.path(input_path, required[1L]), stringsAsFactors = FALSE)
parameters$model_id <- "benchmark"
parameters$model_label <- "Benchmark"
write.csv(parameters, file.path(output_path, "policy-model-parameters.csv"), row.names = FALSE)
saveRDS(
  list(parameters = parameters, cluster_shock = NULL, cluster_external_id = NULL),
  file.path(output_path, "policy-model-parameters.rds"), compress = FALSE
)
write.csv(data.frame(
  model_id = "benchmark", model_label = "Benchmark", model_family = "gaussian",
  lambda_structure = "common", lambda_prior_sd = 0.25,
  draws = nrow(parameters), chains = length(unique(parameters$chain)),
  cluster_shocks = 0L
), file.path(output_path, "policy-model-parameter-status.csv"), row.names = FALSE)
message("Standardized ", nrow(parameters), " balanced benchmark draws.")
