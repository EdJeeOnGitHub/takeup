#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(purrr)
  library(rlang)
})

workspace_path <- option_value(
  "--workspace", "data/stan_analysis_data/dist_fit105.RData"
)
model_name <- option_value(
  "--model", "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)
stan_path <- option_value("--stan-path", "stan_models")
stan_file <- option_value(
  "--stan-file", "takeup_struct_main_core_compact_gq.stan"
)
fit_csvs <- strsplit(
  option_value("--fit-csvs", ""), ",", fixed = TRUE
)[[1L]]
output_path <- option_value("--output-path", "temp/main-core-compact-gq")
use_cluster_shock <- as.integer(
  option_value("--use-core-cluster-shock", "0")
)
shock_prior <- as.numeric(
  option_value("--core-cluster-shock-sd-prior", "0.1")
)
cluster_weight_file <- option_value("--cluster-weight-file")
threads <- as.integer(option_value("--threads", "2"))
parallel_chains <- as.integer(option_value("--parallel-chains", "4"))
cmdstan_path_option <- option_value(
  "--cmdstan-path", Sys.getenv("CMDSTAN", unset = "")
)

if (length(fit_csvs) < 1L || any(!nzchar(fit_csvs)) ||
    any(!file.exists(fit_csvs))) {
  stop("--fit-csvs must list existing fitted-parameter CSVs.", call. = FALSE)
}
if (nzchar(cmdstan_path_option)) {
  Sys.setenv(CMDSTAN = cmdstan_path_option)
  set_cmdstan_path(cmdstan_path_option)
}
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE")
options(cmdstanr_no_ver_check = TRUE)

data <- prepare_main_core_data(
  workspace_path,
  model_name,
  weight_file = cluster_weight_file,
  use_cluster_shock = use_cluster_shock,
  cluster_shock_sd_prior = shock_prior
)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
model <- cmdstan_model(
  file.path(stan_path, stan_file),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
gq_fit <- model$generate_quantities(
  fitted_params = fit_csvs,
  data = data,
  output_dir = output_path,
  output_basename = if (use_cluster_shock) {
    paste0("main-core-shock-", format(shock_prior, trim = TRUE), "-compact")
  } else {
    "main-core-baseline-compact"
  },
  parallel_chains = min(parallel_chains, length(fit_csvs)),
  threads_per_chain = threads
)
message(paste(gq_fit$output_files(), collapse = "\n"))
