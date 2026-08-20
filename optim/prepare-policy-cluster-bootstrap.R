#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")

weighted_path <- policy_option_value(
  args, "--weighted-path",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-exponential-cluster-weights"
)
output_path <- policy_option_value(
  args, "--output-path",
  "optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots-exponential-cluster-weights"
)
num_replicates <- as.integer(policy_option_value(args, "--num-replicates", "999"))
method <- policy_option_value(args, "--method", "exponential")
if (!method %in% c("exponential", "multinomial")) {
  stop("--method must be exponential or multinomial.", call. = FALSE)
}

status_files <- Sys.glob(file.path(weighted_path, "*", "status.csv"))
if (length(status_files) == 0L) stop("No bootstrap status files found.", call. = FALSE)
statuses <- do.call(rbind, lapply(status_files, read.csv, stringsAsFactors = FALSE))
statuses <- statuses[
  statuses$method == method & statuses$status == "complete" &
    file.exists(statuses$mode_csv),
]
statuses <- statuses[order(statuses$replicate), ]
statuses <- statuses[!duplicated(statuses$replicate), ]
parameter_rows <- list()
selected_indices <- integer()
excluded <- data.frame(replicate = integer(), mode_csv = character(), reason = character())
for (index in seq_len(nrow(statuses))) {
  result <- tryCatch(
    canonical_policy_parameters(
      read_cmdstan_mode_row(statuses$mode_csv[index]),
      replicate = statuses$replicate[index],
      mode_csv = statuses$mode_csv[index]
    ),
    error = identity
  )
  if (inherits(result, "error")) {
    excluded <- rbind(excluded, data.frame(
      replicate = statuses$replicate[index],
      mode_csv = statuses$mode_csv[index],
      reason = conditionMessage(result),
      stringsAsFactors = FALSE
    ))
    next
  }
  parameter_rows[[length(parameter_rows) + 1L]] <- result
  selected_indices <- c(selected_indices, index)
  if (length(parameter_rows) == num_replicates) break
}
if (length(parameter_rows) < num_replicates) {
  stop("Only ", length(parameter_rows), " valid bootstrap modes found.", call. = FALSE)
}
parameters <- do.call(rbind, parameter_rows)
statuses <- statuses[selected_indices, ]
parameters$draw <- seq_len(nrow(parameters))

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(parameters, file.path(output_path, "policy-bootstrap-parameters.csv"), row.names = FALSE)
write.csv(excluded, file.path(output_path, "policy-bootstrap-excluded-modes.csv"), row.names = FALSE)
write.csv(
  statuses[, intersect(c(
    "replicate", "seed", "attempts", "elapsed_seconds", "weight_file",
    "mode_csv", "gq_csv"
  ), names(statuses)), drop = FALSE],
  file.path(output_path, "policy-bootstrap-draw-map.csv"), row.names = FALSE
)
message("Prepared ", nrow(parameters), " ", method,
        " cluster-weighted policy modes.")
