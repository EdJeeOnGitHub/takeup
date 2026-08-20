#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

input_path <- option_value(
  "--input-path",
  "/project/akaring/takeup-data/temp-data/indiv-fp-reduce-sum-production-gq"
)
output_path <- option_value("--output-path", input_path)

csv_files <- list.files(
  input_path,
  pattern = "^compact_sm_fit105_.*[.]csv$",
  full.names = TRUE
)
if (length(csv_files) < 2L) {
  stop("Expected at least two compact GQ CSVs in ", input_path, call. = FALSE)
}

# Compact GQ files contain no sampler columns, so reading them directly avoids
# loading CmdStanR and its dependency stack for this lightweight summary job.
draws <- do.call(rbind, lapply(csv_files, function(path) {
  read.csv(path, comment.char = "#", check.names = FALSE)
}))

draw_value <- function(stem, row, column) {
  candidates <- c(
    sprintf("%s[%d,%d]", stem, row, column),
    sprintf("%s.%d.%d", stem, row, column)
  )
  name <- candidates[candidates %in% names(draws)][1L]
  if (is.na(name)) stop("Missing generated quantity: ", candidates[[1L]], call. = FALSE)
  draws[[name]]
}

summarize_value <- function(value) {
  value <- value[is.finite(value)]
  data.frame(
    n = length(value),
    mean = mean(value),
    median = median(value),
    q025 = unname(quantile(value, 0.025)),
    q975 = unname(quantile(value, 0.975))
  )
}

treatments <- c("Control", "Ink", "Calendar", "Bracelet")
levels <- c("Combined", "Close", "Far")
level_draw <- function(level, treatment) {
  draw_value("compact_takeup_level", level, treatment)
}

ate_rows <- list()
for (treatment_index in seq_along(treatments)) {
  for (level_index in seq_along(levels)) {
    value <- if (treatment_index == 1L) {
      level_draw(level_index, 1L)
    } else {
      level_draw(level_index, treatment_index) - level_draw(level_index, 1L)
    }
    ate_rows[[length(ate_rows) + 1L]] <- cbind(
      data.frame(
        estimand = if (treatment_index == 1L) "Control mean" else "ATE vs Control",
        treatment = treatments[[treatment_index]],
        distance_group = levels[[level_index]]
      ),
      summarize_value(value)
    )
  }

  far_close <- if (treatment_index == 1L) {
    level_draw(3L, 1L) - level_draw(2L, 1L)
  } else {
    (level_draw(3L, treatment_index) - level_draw(3L, 1L)) -
      (level_draw(2L, treatment_index) - level_draw(2L, 1L))
  }
  ate_rows[[length(ate_rows) + 1L]] <- cbind(
    data.frame(
      estimand = if (treatment_index == 1L) "Control Far - Close" else "ATE Far - Close",
      treatment = treatments[[treatment_index]],
      distance_group = "Far - Close"
    ),
    summarize_value(far_close)
  )
}
ate_summary <- do.call(rbind, ate_rows)

sm_rows <- list()
distances <- c(500L, 1500L, 2500L)
for (distance_index in seq_along(distances)) {
  for (treatment_index in seq_along(treatments)) {
    # Paper convention: sign-reversed, distance-normalized multiplier.
    value <- -draw_value("compact_sm_rescaled", distance_index, treatment_index)
    sm_rows[[length(sm_rows) + 1L]] <- cbind(
      data.frame(
        treatment = treatments[[treatment_index]],
        distance_m = distances[[distance_index]]
      ),
      summarize_value(value)
    )
  }
}
sm_summary <- do.call(rbind, sm_rows)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(ate_summary, file.path(output_path, "two-chain-ate-summary.csv"), row.names = FALSE)
write.csv(
  sm_summary,
  file.path(output_path, "two-chain-social-multiplier-summary.csv"),
  row.names = FALSE
)

print(ate_summary, row.names = FALSE)
print(sm_summary, row.names = FALSE)
