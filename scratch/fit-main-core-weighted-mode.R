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
  library(posterior)
})

workspace_path <- option_value(
  "--workspace", "data/stan_analysis_data/dist_fit105.RData"
)
data_json <- option_value("--data-json")
model_name <- option_value(
  "--model", "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)
stan_path <- option_value("--stan-path", "stan_models")
stan_file <- option_value("--stan-file", "takeup_struct_main_core.stan")
gq_file <- option_value(
  "--gq-file", "takeup_struct_main_core_compact_gq.stan"
)
output_path <- option_value("--output-path", "temp/main-core-weighted-modes")
label_override <- option_value("--label")
weight_file <- option_value("--weight-file")
manifest_path <- option_value("--manifest")
task_id <- as.integer(option_value("--task-id", "0"))
init_json <- option_value("--init-json")
threads <- as.integer(option_value("--threads", "8"))
iterations <- as.integer(option_value("--iterations", "2000"))
seed <- as.integer(option_value("--seed", "20260802"))
cmdstan_path_option <- option_value(
  "--cmdstan-path", Sys.getenv("CMDSTAN", unset = "")
)
lambda_structure <- as.integer(option_value("--core-lambda-structure", "0"))
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

method <- "unweighted"
replicate_id <- 0L
if (!is.null(manifest_path)) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  if (task_id < 1L || task_id > nrow(manifest)) {
    stop("--task-id is outside the manifest.", call. = FALSE)
  }
  task <- manifest[task_id, ]
  weight_file <- task$weight_file
  method <- task$method
  replicate_id <- as.integer(task$replicate)
  seed <- seed + task_id
}
if (!is.null(weight_file) && method == "unweighted") method <- "weighted"
label <- if (!is.null(label_override)) {
  label_override
} else if (replicate_id > 0L) {
  sprintf("%s-%04d", method, replicate_id)
} else {
  method
}

if (!is.null(cmdstan_path_option) && nzchar(cmdstan_path_option)) {
  Sys.setenv(CMDSTAN = cmdstan_path_option)
  set_cmdstan_path(cmdstan_path_option)
}
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE", STAN_NUM_THREADS = threads)
options(cmdstanr_no_ver_check = TRUE)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
run_path <- file.path(output_path, label)
dir.create(run_path, recursive = TRUE, showWarnings = FALSE)

data <- if (!is.null(data_json)) {
  if (!file.exists(data_json)) {
    stop("Data JSON not found: ", data_json, call. = FALSE)
  }
  fitted_data_json <- file.path(run_path, "fit-data.json")
  main_core_patch_stan_json_scalars(
    data_json,
    fitted_data_json,
    c(
      core_lambda_structure = lambda_structure,
      core_lambda_log_ratio_sd_prior = lambda_prior,
      core_profile_group_lambda = profile_group_lambda,
      core_profile_group_log_ratio = profile_group_log_ratio,
      core_type_distribution = type_distribution,
      core_student_t_df = student_t_df,
      core_student_t_components = student_t_components
    )
  )
} else {
  prepare_main_core_data(
    workspace_path = workspace_path,
    model_name = model_name,
    weight_file = weight_file,
    use_cluster_shock = 0L,
    lambda_structure = lambda_structure,
    lambda_log_ratio_sd_prior = lambda_prior,
    profile_group_lambda = profile_group_lambda,
    profile_group_log_ratio = profile_group_log_ratio,
    type_distribution = type_distribution,
    student_t_df = student_t_df,
    student_t_components = student_t_components
  )
}
model <- cmdstan_model(
  file.path(stan_path, stan_file),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
gq_model <- cmdstan_model(
  file.path(stan_path, gq_file),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)

init <- if (is.null(init_json)) 0 else init_json
started <- Sys.time()
attempt <- 1L
error_message <- NA_character_
fit <- tryCatch(
  model$optimize(
    data = data,
    seed = seed,
    init = init,
    algorithm = "lbfgs",
    iter = iterations,
    threads = threads,
    output_dir = run_path,
    output_basename = paste0("mode-", label),
    refresh = 100
  ),
  error = identity
)
if (inherits(fit, "error")) {
  attempt <- 2L
  error_message <- conditionMessage(fit)
  fit <- tryCatch(
    model$optimize(
      data = data,
      seed = seed,
      init = init,
      algorithm = "lbfgs",
      iter = 2L * iterations,
      init_alpha = 0.01,
      threads = threads,
      output_dir = run_path,
      output_basename = paste0("mode-", label, "-retry"),
      refresh = 100
    ),
    error = identity
  )
}
if (inherits(fit, "error")) {
  attempt <- 3L
  error_message <- paste(
    na.omit(c(error_message, conditionMessage(fit))), collapse = " | "
  )
  fit <- tryCatch(
    model$optimize(
      data = data,
      seed = seed,
      init = init,
      algorithm = "newton",
      iter = 200L,
      threads = threads,
      output_dir = run_path,
      output_basename = paste0("mode-", label, "-newton"),
      refresh = 10
    ),
    error = identity
  )
}

status <- if (inherits(fit, "error")) "failed" else "complete"
if (status == "complete") {
  mode_csv <- fit$output_files()
  init_output <- file.path(run_path, "mode-init.json")
  write_mode_init_json(fit, model, init_output)
  gq_fit <- gq_model$generate_quantities(
    fitted_params = mode_csv,
    data = data,
    seed = seed,
    output_dir = run_path,
    output_basename = paste0("compact-gq-", label),
    parallel_chains = 1,
    threads_per_chain = threads
  )
  gq_csv <- gq_fit$output_files()
} else {
  mode_csv <- NA_character_
  init_output <- NA_character_
  gq_csv <- NA_character_
  error_message <- paste(na.omit(c(error_message, conditionMessage(fit))), collapse = " | ")
}

status_data <- data.frame(
  method = method,
  label = label,
  lambda_structure = lambda_structure,
  lambda_log_ratio_sd_prior = lambda_prior,
  profile_group_lambda = profile_group_lambda,
  profile_group_log_ratio = profile_group_log_ratio,
  type_distribution = type_distribution,
  student_t_df = student_t_df,
  student_t_components = student_t_components,
  replicate = replicate_id,
  seed = seed,
  status = status,
  attempts = attempt,
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  weight_file = if (is.null(weight_file)) NA_character_ else weight_file,
  mode_csv = mode_csv,
  init_json = init_output,
  gq_csv = gq_csv,
  error = error_message,
  stringsAsFactors = FALSE
)
write.csv(status_data, file.path(run_path, "status.csv"), row.names = FALSE)
if (status != "complete") stop(error_message, call. = FALSE)
print(status_data)
