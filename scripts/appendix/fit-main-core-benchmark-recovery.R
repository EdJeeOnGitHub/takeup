#!/usr/bin/env Rscript

# Fit one benchmark-recovery dataset. A failed diagnostic audit automatically
# triggers a longer HMC attempt in the same job.

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

manifest_path <- option_value("--manifest")
task_id <- as.integer(option_value("--task-id", "0"))
output_path <- option_value(
  "--output-path", "temp-data/main-core-benchmark-recovery"
)
stan_path <- option_value("--stan-path", "stan_models")
cmdstan_path_option <- option_value(
  "--cmdstan-path", Sys.getenv("CMDSTAN_PATH", unset = "")
)
threads <- as.integer(option_value("--threads", "8"))
chains <- as.integer(option_value("--chains", "1"))
iter_warmup <- as.integer(option_value("--iter-warmup", "1000"))
iter_sampling <- as.integer(option_value("--iter-sampling", "1000"))
rerun_warmup <- as.integer(option_value("--rerun-warmup", "2000"))
rerun_sampling <- as.integer(option_value("--rerun-sampling", "2000"))
minimum_ess <- as.numeric(option_value("--minimum-ess", "400"))
maximum_rhat <- as.numeric(option_value("--maximum-rhat", "1.01"))
enforce_minimum_ess <- tolower(option_value(
  "--enforce-minimum-ess", "true"
)) %in% c("1", "true", "yes")
if (chains < 1L || chains > threads) stop("--chains must be between 1 and --threads.")

if (is.null(manifest_path) || !file.exists(manifest_path)) {
  stop("--manifest must identify the recovery manifest.")
}
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(manifest)) stop("Invalid --task-id.")
task <- manifest[task_id, ]
truth_csv <- file.path(output_path, "truth-parameters.csv")
candidate_csv <- file.path(output_path, "likelihood-candidates.csv")
if (!all(file.exists(c(task$data_json, truth_csv, candidate_csv)))) {
  stop("Recovery inputs are incomplete for task ", task_id)
}
if (nzchar(cmdstan_path_option)) set_cmdstan_path(cmdstan_path_option)
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE", STAN_NUM_THREADS = threads)
options(cmdstanr_no_ver_check = TRUE)

model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core.stan"),
  include_paths = stan_path, cpp_options = list(stan_threads = TRUE)
)
gq_model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core_compact_gq.stan"),
  include_paths = stan_path, cpp_options = list(stan_threads = TRUE)
)
data <- jsonlite::read_json(task$data_json, simplifyVector = TRUE)
if (length(data$optim_distances) == 1L) {
  data$optim_distances <- array(data$optim_distances, dim = 1L)
}
gq_data <- data
gq_data$roc_distances[[6L]] <- gq_data$roc_distances[[6L]] * 2.5

reported_parameters <- c(
  "dist_beta_v[1]", "base_mu_rep", "raw_u_sd[1]", "wtp_value_utility"
)
reported_multipliers <- paste0("core_compact_sm_rescaled[1,", 1:4, "]")
run_root <- file.path(output_path, "fits", sprintf("rep-%03d", task$replicate))
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

run_attempt <- function(attempt, warmup, sampling) {
  attempt_path <- file.path(run_root, paste0("attempt-", attempt))
  dir.create(attempt_path, recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()
  mode <- model$optimize(
    data = data, seed = task$seed + attempt * 1000L, init = 0,
    algorithm = "lbfgs", iter = 2000, threads = threads,
    output_dir = attempt_path, output_basename = "posterior-mode", refresh = 100
  )
  init_json <- file.path(attempt_path, "mode-init.json")
  write_mode_init_json(mode, model, init_json)
  fit <- model$sample(
    data = data, seed = task$seed + attempt * 10000L,
    chains = chains, parallel_chains = chains,
    threads_per_chain = max(1L, threads %/% chains),
    iter_warmup = warmup, iter_sampling = sampling, init = init_json,
    metric = "diag_e", adapt_delta = 0.999, max_treedepth = 12,
    output_dir = attempt_path, output_basename = "benchmark-recovery",
    refresh = 100
  )
  gq <- gq_model$generate_quantities(
    fitted_params = fit$output_files(), data = gq_data,
    seed = task$seed + attempt * 20000L, output_dir = attempt_path,
    output_basename = "compact-gq-1250m", parallel_chains = chains,
    threads_per_chain = max(1L, threads %/% chains)
  )
  parameter_summary <- fit$summary(variables = reported_parameters)
  gq_summary <- summarise_draws(
    read_cmdstan_csv(gq$output_files(), variables = reported_multipliers)$generated_quantities
  )
  sampler <- as_draws_df(fit$sampler_diagnostics())
  max_rhat <- if (chains == 1L) NA_real_ else
    max(c(parameter_summary$rhat, gq_summary$rhat), na.rm = TRUE)
  audit <- data.frame(
    attempt = attempt, warmup = warmup, sampling = sampling,
    max_rhat = max_rhat,
    min_ess_bulk = min(c(parameter_summary$ess_bulk, gq_summary$ess_bulk), na.rm = TRUE),
    min_ess_tail = min(c(parameter_summary$ess_tail, gq_summary$ess_tail), na.rm = TRUE),
    divergences = sum(sampler$divergent__, na.rm = TRUE),
    max_treedepth = sum(sampler$treedepth__ >= 12, na.rm = TRUE),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  audit$ess_advisory_passed <- with(audit,
    min_ess_bulk >= minimum_ess && min_ess_tail >= minimum_ess
  )
  audit$ess_enforced <- enforce_minimum_ess
  audit$passed <- with(audit,
    (chains == 1L || (is.finite(max_rhat) && max_rhat <= maximum_rhat)) &&
      (!ess_enforced || ess_advisory_passed) &&
      divergences == 0 && max_treedepth == 0
  )
  write.csv(audit, file.path(attempt_path, "audit.csv"), row.names = FALSE)
  list(fit = fit, gq = gq, audit = audit, path = attempt_path)
}

result <- tryCatch(
  run_attempt(1L, iter_warmup, iter_sampling),
  error = function(error) structure(list(error = conditionMessage(error)), class = "failed_attempt")
)
attempts <- list(result)
if (inherits(result, "failed_attempt") || !isTRUE(result$audit$passed)) {
  result <- tryCatch(
    run_attempt(2L, rerun_warmup, rerun_sampling),
    error = function(error) structure(list(error = conditionMessage(error)), class = "failed_attempt")
  )
  attempts[[2L]] <- result
}

passed <- !inherits(result, "failed_attempt") && isTRUE(result$audit$passed)
likelihood_file <- NA_character_
if (passed) {
  likelihood_gq <- gq_model$generate_quantities(
    fitted_params = candidate_csv, data = gq_data, seed = task$seed + 30000L,
    output_dir = result$path, output_basename = "likelihood-slices",
    threads_per_chain = threads
  )
  likelihood_file <- likelihood_gq$output_files()[1]
}
status <- data.frame(
  task_id = task_id, replicate = task$replicate, chains = chains,
  status = if (passed) "complete" else "failed",
  final_attempt = length(attempts),
  fit_csvs = if (passed) paste(result$fit$output_files(), collapse = ";") else NA_character_,
  gq_csvs = if (passed) paste(result$gq$output_files(), collapse = ";") else NA_character_,
  likelihood_gq = likelihood_file,
  error = if (inherits(result, "failed_attempt")) result$error else NA_character_,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write.csv(status, file.path(run_root, "status.csv"), row.names = FALSE)
if (!passed) stop("Recovery task ", task_id, " failed its HMC audit.")
message("Recovery task ", task_id, " complete after ", length(attempts), " attempt(s).")
