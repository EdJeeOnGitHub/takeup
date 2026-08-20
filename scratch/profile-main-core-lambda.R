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

manifest_path <- option_value("--manifest")
task_id <- as.integer(option_value("--task-id", "0"))
workspace <- option_value("--workspace")
model_name <- option_value(
  "--model", "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)
stan_path <- option_value("--stan-path", "stan_models")
output_path <- option_value(
  "--output-path", "temp-data/main-core-lambda-identification/profile"
)
init_json <- option_value("--init-json")
threads <- as.integer(option_value("--threads", "8"))
cmdstan_path_option <- option_value(
  "--cmdstan-path", Sys.getenv("CMDSTAN", unset = "")
)
if (is.null(manifest_path) || !file.exists(manifest_path) ||
    is.null(workspace) || !file.exists(workspace)) {
  stop("Existing --manifest and --workspace files are required.", call. = FALSE)
}
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(manifest)) {
  stop("--task-id is outside the profile manifest.", call. = FALSE)
}
task <- manifest[task_id, ]
if (nzchar(cmdstan_path_option)) {
  Sys.setenv(CMDSTAN = cmdstan_path_option)
  set_cmdstan_path(cmdstan_path_option)
}
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE", STAN_NUM_THREADS = threads)
options(cmdstanr_no_ver_check = TRUE)

data <- prepare_main_core_data(
  workspace,
  model_name,
  lambda_structure = 1L,
  lambda_log_ratio_sd_prior = 0.25,
  profile_group_lambda = 1L,
  profile_group_log_ratio = task$theta
)
model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core.stan"),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
gq_model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core_compact_gq.stan"),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
run_path <- file.path(output_path, "fits", task$label)
dir.create(run_path, recursive = TRUE, showWarnings = FALSE)

profile_init <- 0
if (!is.null(init_json) && file.exists(init_json)) {
  init_values <- jsonlite::read_json(init_json, simplifyVector = TRUE)
  init_values$core_lambda_group_log_ratio_raw <- NULL
  profile_init <- file.path(run_path, "profile-init.json")
  write_stan_json(init_values, profile_init)
}

initializers <- list(
  mode = profile_init,
  zero = 0,
  random_small = 0.05
)
fits <- vector("list", length(initializers))
attempt_status <- vector("list", length(initializers))
for (index in seq_along(initializers)) {
  attempt_name <- names(initializers)[index]
  started <- Sys.time()
  fit <- tryCatch(
    model$optimize(
      data = data,
      init = initializers[[index]],
      seed = task$seed + index,
      algorithm = "lbfgs",
      iter = 3000,
      threads = threads,
      output_dir = run_path,
      output_basename = paste0("profile-", task$label, "-", attempt_name),
      refresh = 250
    ),
    error = identity
  )
  lp_result <- if (inherits(fit, "error")) fit else tryCatch(
    as.numeric(as_draws_df(fit$draws(variables = "lp__"))$lp__)[1],
    error = identity
  )
  lp <- if (inherits(lp_result, "error")) NA_real_ else lp_result
  fits[[index]] <- fit
  attempt_status[[index]] <- data.frame(
    task_id = task_id,
    theta = task$theta,
    attempt = attempt_name,
    success = !inherits(fit, "error") && is.finite(lp),
    lp = lp,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    error = if (inherits(lp_result, "error")) {
      conditionMessage(lp_result)
    } else {
      NA_character_
    }
  )
}
attempts <- bind_rows(attempt_status)
write.csv(attempts, file.path(run_path, "attempts.csv"), row.names = FALSE)
if (all(!attempts$success) || all(!is.finite(attempts$lp))) {
  stop("All profile optimizations failed at theta=", task$theta, call. = FALSE)
}
best_index <- which.max(ifelse(is.finite(attempts$lp), attempts$lp, -Inf))
best_fit <- fits[[best_index]]
best_csv <- best_fit$output_files()
gq <- gq_model$generate_quantities(
  fitted_params = best_csv,
  data = data,
  seed = task$seed,
  output_dir = run_path,
  output_basename = paste0("profile-gq-", task$label),
  threads_per_chain = threads
)
status <- data.frame(
  task_id = task_id,
  theta = task$theta,
  label = task$label,
  status = "complete",
  best_attempt = attempts$attempt[best_index],
  lp = attempts$lp[best_index],
  fit_csv = best_csv,
  gq_csv = gq$output_files()
)
write.csv(status, file.path(run_path, "status.csv"), row.names = FALSE)
print(status)
