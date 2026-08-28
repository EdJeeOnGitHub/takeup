#!/usr/bin/env Rscript

# Validate the policy postprocessing rebuilt from streamlined active-robustness
# fits. This is a postprocessing audit only; it never samples a Stan model.

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit)) substring(hit, nchar(prefix) + 1L) else default
}

optim_root <- "/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
new_root <- value(
  "--new-root", file.path(optim_root, "policy-model-robustness-streamlined")
)
old_root <- value(
  "--old-root", file.path(optim_root, "policy-model-robustness")
)
output_path <- value(
  "--output-path", file.path(new_root, "policy-migration-audit.csv")
)
comparison_path <- value(
  "--comparison-path", file.path(new_root, "policy-summary-comparison.csv")
)
expected_draws <- as.integer(value("--expected-draws", "4000"))

models <- c("exclude-dispersed", "second-order-observability")
scenarios <- c(
  "control", "bracelet", "static-control", "static-bracelet",
  "suppress-reputation"
)

read_one <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

audit_model <- function(model_id) {
  directory <- file.path(new_root, model_id)
  required <- c(
    "policy-model-summary.csv", "policy-model-contrast-summary.csv",
    "policy-model-parameter-status.csv", "policy-prediction-status.csv",
    "policy-model-parameters.csv", "provenance.txt"
  )
  missing <- required[!file.exists(file.path(directory, required))]
  parameter <- if (!length(missing)) {
    read_one(file.path(directory, "policy-model-parameter-status.csv"))
  } else data.frame()
  prediction <- if (!length(missing)) {
    read_one(file.path(directory, "policy-prediction-status.csv"))
  } else data.frame()
  summary <- if (!length(missing)) {
    read_one(file.path(directory, "policy-model-summary.csv"))
  } else data.frame()

  status_rows <- integer(length(scenarios))
  status_draws <- integer(length(scenarios))
  allocation_files <- integer(length(scenarios))
  error_rows <- integer(length(scenarios))
  for (index in seq_along(scenarios)) {
    scenario <- scenarios[[index]]
    allocation_dir <- file.path(directory, "allocations", scenario)
    status_path <- file.path(allocation_dir, "status.csv")
    if (!file.exists(status_path)) next
    status <- read.csv(status_path, stringsAsFactors = FALSE)
    status_rows[[index]] <- nrow(status)
    status_draws[[index]] <- length(unique(status$draw))
    allocation_files[[index]] <- length(list.files(
      allocation_dir, pattern = "^replicate-[0-9]{4}[.]rds$"
    ))
    error_rows[[index]] <- if ("status" %in% names(status)) {
      sum(status$status %in% c("error", "failed"))
    } else nrow(status)
  }

  provenance_path <- file.path(directory, "provenance.txt")
  provenance <- if (file.exists(provenance_path)) {
    readLines(provenance_path, warn = FALSE)
  } else character()
  streamlined_provenance <- any(
    provenance == "streamlined_active_robustness=true"
  )
  correct_belief_order <- model_id != "second-order-observability" || any(
    provenance == "extract_options=--beliefs-order 2"
  )
  undefined <- if (nrow(prediction) &&
      "undefined_prediction_values" %in% names(prediction)) {
    prediction$undefined_prediction_values[[1L]]
  } else NA_real_

  pass <-
    !length(missing) &&
    nrow(parameter) == 1L && parameter$draws[[1L]] == expected_draws &&
    parameter$chains[[1L]] == 4L &&
    nrow(prediction) == 1L && prediction$draws[[1L]] == expected_draws &&
    prediction$scenarios[[1L]] == length(scenarios) &&
    is.finite(undefined) && undefined == 0 &&
    nrow(summary) == length(scenarios) + 1L &&
    all(summary$draws == expected_draws) &&
    all(status_rows == expected_draws) &&
    all(status_draws == expected_draws) &&
    all(allocation_files == expected_draws) &&
    sum(error_rows) == 0L &&
    streamlined_provenance && correct_belief_order

  data.frame(
    model_id = model_id,
    expected_draws = expected_draws,
    parameter_draws = if (nrow(parameter)) parameter$draws[[1L]] else NA,
    chains = if (nrow(parameter)) parameter$chains[[1L]] else NA,
    prediction_draws = if (nrow(prediction)) prediction$draws[[1L]] else NA,
    undefined_prediction_values = undefined,
    summary_rows = nrow(summary),
    minimum_scenario_status_rows = min(status_rows),
    minimum_scenario_unique_draws = min(status_draws),
    minimum_scenario_allocation_files = min(allocation_files),
    error_status_rows = sum(error_rows),
    streamlined_provenance = streamlined_provenance,
    correct_belief_order = correct_belief_order,
    missing_files = paste(missing, collapse = ";"),
    pass = pass,
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, lapply(models, audit_model))
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(audit, output_path, row.names = FALSE, quote = TRUE)

comparisons <- lapply(models, function(model_id) {
  new_path <- file.path(new_root, model_id, "policy-model-summary.csv")
  old_path <- file.path(old_root, model_id, "policy-model-summary.csv")
  if (!file.exists(new_path) || !file.exists(old_path)) return(NULL)
  new <- read.csv(new_path, stringsAsFactors = FALSE)
  old <- read.csv(old_path, stringsAsFactors = FALSE)
  columns <- intersect(
    c("n_pot_estimate", "takeup_estimate", "distance_estimate"),
    intersect(names(new), names(old))
  )
  merged <- merge(
    old[c("scenario", columns)], new[c("scenario", columns)],
    by = "scenario", suffixes = c("_old", "_streamlined")
  )
  for (column in columns) {
    merged[[paste0(column, "_difference")]] <-
      merged[[paste0(column, "_streamlined")]] -
      merged[[paste0(column, "_old")]]
  }
  data.frame(model_id = model_id, merged, stringsAsFactors = FALSE)
})
comparisons <- Filter(Negate(is.null), comparisons)
comparison <- if (length(comparisons)) do.call(rbind, comparisons) else data.frame()
write.csv(comparison, comparison_path, row.names = FALSE, quote = TRUE)

if (!all(audit$pass)) {
  stop("At least one streamlined policy replacement failed its audit.")
}
message("All streamlined policy replacement audits passed: ", output_path)
