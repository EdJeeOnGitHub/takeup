#!/usr/bin/env Rscript

args <- docopt::docopt(
  "
Generate blockwise log-likelihood totals for an existing structural fit.

Usage:
  scripts/appendix/generate_bayesian_sensitivity_quantities.R <fit-version> --model=<model> [options]

Options:
  --analysis-file=<path>  Saved scripts/structural/run-model.R workspace [default: data/stan_analysis_data/dist_fit{fit-version}.RData]
  --input-path=<path>     Directory containing fitted CmdStan CSV files [default: data/stan_analysis_data]
  --output-path=<path>    Directory for generated-quantities CSV files [default: temp-data/bayesian-sensitivity]
  --include-path=<path>   Stan include directory [default: stan_models]
  --chains=<chains>       Comma-separated chain numbers [default: 1,2,3,4]
  --seed=<n>              Generated-quantities seed [default: 8675309]
  ",
  args = commandArgs(trailingOnly = TRUE)
)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(purrr)
  library(tidyverse)
})

source("R/structural/legacy-utils.R")

interpolate_fit_version <- function(path, fit_version) {
  stringr::str_replace_all(path, fixed("{fit-version}"), fit_version)
}

fit_version <- args$fit_version
model_name <- args$model
analysis_file <- interpolate_fit_version(args$analysis_file, fit_version)
chains <- as.integer(strsplit(args$chains, ",", fixed = TRUE)[[1]])

if (!file.exists(analysis_file)) {
  stop("Analysis workspace does not exist: ", analysis_file)
}

fit_files <- file.path(
  args$input_path,
  sprintf("dist_fit%s_%s-%d.csv", fit_version, model_name, chains)
)

if (any(!file.exists(fit_files))) {
  stop(
    "Missing fitted CmdStan CSV files:\n",
    paste(fit_files[!file.exists(fit_files)], collapse = "\n")
  )
}

fit_environment <- new.env(parent = globalenv())
load(analysis_file, envir = fit_environment)

if (!exists("models", envir = fit_environment, inherits = FALSE) ||
    !exists("stan_data", envir = fit_environment, inherits = FALSE)) {
  stop("The analysis workspace must contain objects named 'models' and 'stan_data'.")
}

models <- get("models", envir = fit_environment)
stan_data <- get("stan_data", envir = fit_environment)

if (!model_name %in% names(models)) {
  stop(
    "Model '", model_name, "' is not present in ", analysis_file,
    ". Available models: ", paste(names(models), collapse = ", ")
  )
}

model_info <- models[[model_name]]
stan_data_preprocess <- model_info$stan_data_preprocess %||% identity
model_info_for_stan <- model_info
model_info_for_stan$stan_data_preprocess <- NULL
model_info_for_stan$model_file <- "takeup_struct_sensitivity.stan"

sensitivity_data <- stan_data_preprocess(stan_data) %>%
  map_at(
    c("cluster_treatment_map", "beliefs_ate_pairs"),
    ~ mutate(.x, across(everything(), as.integer)) %>% as.matrix()
  ) %>%
  list_modify(!!!model_info_for_stan) %>%
  map_if(is.factor, as.integer) %>%
  discard(~ is_function(.x) || is.character(.x))

sensitivity_data$analysis_data <- NULL

# Fit 104 predates the explicit missing-beliefs indicators. In that fit every
# row entering the beliefs submodel was observed, and the cluster-row shortcut
# was disabled. Supply the equivalent values required by the current Stan data
# declaration so that old posterior draws can be processed reproducibly.
if (is.null(sensitivity_data$use_belief_row_cluster_mu_rep)) {
  sensitivity_data$use_belief_row_cluster_mu_rep <- 0L
}
if (is.null(sensitivity_data$num_belief_rows_by_cluster)) {
  sensitivity_data$num_belief_rows_by_cluster <-
    integer(sensitivity_data$num_clusters)
}
if (is.null(sensitivity_data$belief_observed)) {
  sensitivity_data$belief_observed <-
    rep.int(1L, sensitivity_data$num_beliefs_obs)
}

dir.create(args$output_path, recursive = TRUE, showWarnings = FALSE)

sensitivity_model <- cmdstan_model(
  file.path("stan_models", "takeup_struct_sensitivity.stan"),
  cpp_options = list(stan_threads = TRUE),
  include_paths = args$include_path
)

message("Generating four blockwise likelihood totals from ", length(fit_files), " chains...")
sensitivity_fit <- sensitivity_model$generate_quantities(
  data = sensitivity_data,
  fitted_params = fit_files,
  seed = as.integer(args$seed),
  parallel_chains = min(length(fit_files), 4),
  threads_per_chain = 1
)

output_basename <- sprintf(
  "dist_fit%s_%s_sensitivity",
  fit_version,
  model_name
)
sensitivity_fit$save_output_files(
  dir = args$output_path,
  basename = output_basename,
  timestamp = FALSE,
  random = FALSE
)

message("Saved sensitivity generated quantities to ", args$output_path)
