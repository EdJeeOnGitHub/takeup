#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
output_path <- main_core_option_value(
  args, "--output-path", "temp-data/main-core-lambda-identification/profile"
)
grid_min <- as.numeric(main_core_option_value(args, "--grid-min", "-3"))
grid_max <- as.numeric(main_core_option_value(args, "--grid-max", "3"))
grid_step <- as.numeric(main_core_option_value(args, "--grid-step", "0.1"))
if (!all(is.finite(c(grid_min, grid_max, grid_step))) ||
    grid_min >= grid_max || grid_step <= 0) {
  stop("Invalid profile grid.", call. = FALSE)
}
theta <- seq(grid_min, grid_max, by = grid_step)
manifest <- data.frame(
  task_id = seq_along(theta),
  theta = theta,
  label = paste0(
    "theta-", ifelse(theta < 0, "m", "p"), sprintf("%.2f", abs(theta))
  ),
  seed = 20262000L + seq_along(theta),
  status = "pending"
)
manifest$label <- gsub("\\.", "p", manifest$label)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(
  manifest,
  file.path(output_path, "lambda-profile-manifest.csv"),
  row.names = FALSE
)
print(manifest, row.names = FALSE)
