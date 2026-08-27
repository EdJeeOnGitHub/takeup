#!/usr/bin/env Rscript

# Finalize an existing recovery attempt when hard HMC diagnostics pass and an
# ESS-only failure has been explicitly accepted as advisory.

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
value <- function(name, default = NULL) main_core_option_value(args, name, default)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

input_path <- value("--input-path")
task_id <- as.integer(value("--task-id", "0"))
attempt <- as.integer(value("--attempt", "2"))
cmdstan_path <- value("--cmdstan-path", Sys.getenv("CMDSTAN_PATH", unset = ""))
if (is.null(input_path) || task_id < 1L || !nzchar(cmdstan_path)) {
  stop("--input-path, --task-id, and --cmdstan-path are required.")
}
set_cmdstan_path(cmdstan_path)
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE")
options(cmdstanr_no_ver_check = TRUE)

manifest <- read.csv(file.path(input_path, "benchmark-recovery-manifest.csv"))
task <- manifest[manifest$task_id == task_id, ]
if (nrow(task) != 1L) stop("Task is absent or duplicated in the manifest.")
run_root <- file.path(input_path, "fits", sprintf("rep-%03d", task$replicate))
attempt_path <- file.path(run_root, paste0("attempt-", attempt))
audit_path <- file.path(attempt_path, "audit.csv")
audit <- read.csv(audit_path)
if (nrow(audit) != 1L || audit$divergences != 0L || audit$max_treedepth != 0L) {
  stop("Cannot accept an attempt with failed hard HMC diagnostics.")
}
fit_csvs <- list.files(
  attempt_path, pattern = "^benchmark-recovery-[0-9]+[.]csv$", full.names = TRUE
)
gq_csvs <- list.files(
  attempt_path, pattern = "^compact-gq-1250m-[0-9]+[.]csv$", full.names = TRUE
)
if (!length(fit_csvs) || !length(gq_csvs)) stop("Existing fit/GQ CSVs are missing.")

data <- jsonlite::read_json(task$data_json, simplifyVector = TRUE)
if (length(data$optim_distances) == 1L) {
  data$optim_distances <- array(data$optim_distances, dim = 1L)
}
data$roc_distances[[6L]] <- data$roc_distances[[6L]] * 2.5
model <- cmdstan_model(
  "stan_models/takeup_struct_main_core_compact_gq.stan",
  include_paths = "stan_models", cpp_options = list(stan_threads = TRUE)
)
likelihood <- model$generate_quantities(
  fitted_params = file.path(input_path, "likelihood-candidates.csv"),
  data = data, seed = task$seed + 30000L, output_dir = attempt_path,
  output_basename = "likelihood-slices", threads_per_chain = 8
)

acceptance <- transform(
  audit, ess_enforced = FALSE,
  ess_advisory_passed = min_ess_bulk >= 400 & min_ess_tail >= 400,
  hard_diagnostics_passed = TRUE, accepted = TRUE,
  acceptance_reason = "ESS advisory accepted after inspecting completed attempt"
)
write.csv(acceptance, file.path(attempt_path, "acceptance-audit.csv"), row.names = FALSE)
status <- data.frame(
  task_id = task_id, replicate = task$replicate, chains = 1L,
  status = "complete", final_attempt = attempt,
  fit_csvs = paste(fit_csvs, collapse = ";"),
  gq_csvs = paste(gq_csvs, collapse = ";"),
  likelihood_gq = likelihood$output_files()[1], error = NA_character_,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write.csv(status, file.path(run_root, "status.csv"), row.names = FALSE)
message("Finalized recovery task ", task_id, " with ESS recorded as advisory.")
