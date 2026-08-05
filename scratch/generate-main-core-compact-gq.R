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
data_json <- option_value("--data-json")
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
output_basename <- option_value("--output-basename")
use_cluster_shock <- as.integer(
  option_value("--use-core-cluster-shock", "0")
)
shock_prior <- as.numeric(
  option_value("--core-cluster-shock-sd-prior", "0.1")
)
lambda_structure <- as.integer(
  option_value("--core-lambda-structure", "0")
)
lambda_prior <- as.numeric(
  option_value("--core-lambda-log-ratio-sd-prior", "0.25")
)
profile_group_lambda <- as.integer(
  option_value("--core-profile-group-lambda", "0")
)
profile_group_log_ratio <- as.numeric(
  option_value("--core-profile-group-log-ratio", "0")
)
type_distribution <- as.integer(option_value("--core-type-distribution", "0"))
student_t_df <- as.numeric(option_value("--core-student-t-df", "5"))
student_t_components <- as.integer(
  option_value("--core-student-t-components", "12")
)
observation_model <- as.integer(
  option_value("--core-observation-model", "0")
)
recognition_structure <- as.integer(
  option_value("--core-recognition-structure", "0")
)
report_structure <- as.integer(
  option_value("--core-report-structure", "0")
)
peer_audit_path <- option_value("--peer-audit-path")
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

data <- if (!is.null(data_json)) {
  if (observation_model > 0L) {
    stop("Asymmetric observability requires rebuilding data, not --data-json.",
         call. = FALSE)
  }
  if (!file.exists(data_json)) stop("Data JSON not found: ", data_json, call. = FALSE)
  patched_data_json <- tempfile(fileext = ".json")
  on.exit(unlink(patched_data_json), add = TRUE)
  main_core_patch_stan_json_scalars(
    data_json,
    patched_data_json,
    c(
      core_lambda_structure = lambda_structure,
      core_lambda_log_ratio_sd_prior = lambda_prior,
      core_profile_group_lambda = profile_group_lambda,
      core_profile_group_log_ratio = profile_group_log_ratio,
      core_gq_override_lambda = 0L
    )
  )
} else {
  prepare_main_core_data(
    workspace_path,
    model_name,
    weight_file = cluster_weight_file,
    use_cluster_shock = use_cluster_shock,
    cluster_shock_sd_prior = shock_prior,
    lambda_structure = lambda_structure,
    lambda_log_ratio_sd_prior = lambda_prior,
    profile_group_lambda = profile_group_lambda,
    profile_group_log_ratio = profile_group_log_ratio,
    type_distribution = type_distribution,
    student_t_df = student_t_df,
    student_t_components = student_t_components,
    observation_model = observation_model,
    recognition_structure = recognition_structure,
    report_structure = report_structure,
    peer_audit_path = peer_audit_path
  )
}
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
  output_basename = if (!is.null(output_basename)) {
    output_basename
  } else if (type_distribution == 1L) {
    paste0("main-core-student-t-df", format(student_t_df, trim = TRUE), "-compact")
  } else if (observation_model > 0L) {
    paste0(
      "main-core-observation-", observation_model, "-rec-",
      recognition_structure, "-report-", report_structure, "-compact"
    )
  } else if (lambda_structure > 0) {
    paste0(
      "main-core-lambda-", lambda_structure, "-",
      format(lambda_prior, trim = TRUE),
      "-compact"
    )
  } else if (use_cluster_shock) {
    paste0("main-core-shock-", format(shock_prior, trim = TRUE), "-compact")
  } else {
    "main-core-baseline-compact"
  },
  parallel_chains = min(parallel_chains, length(fit_csvs)),
  threads_per_chain = threads
)
message(paste(gq_fit$output_files(), collapse = "\n"))
