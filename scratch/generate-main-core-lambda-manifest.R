#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
output_path <- main_core_option_value(
  args,
  "--output-path",
  "temp-data/main-core-lambda-identification"
)
prior_grid <- as.numeric(strsplit(
  main_core_option_value(args, "--prior-grid", "0.10,0.25,0.50,1.00"),
  ",",
  fixed = TRUE
)[[1L]])
if (length(prior_grid) < 1L || any(!is.finite(prior_grid)) ||
    any(prior_grid <= 0) || anyDuplicated(prior_grid)) {
  stop("--prior-grid must contain distinct positive numbers.", call. = FALSE)
}

specifications <- rbind(
  data.frame(
    structure = "common", lambda_structure = 0L,
    lambda_log_ratio_sd_prior = 0.25
  ),
  data.frame(
    structure = "grouped", lambda_structure = 1L,
    lambda_log_ratio_sd_prior = prior_grid
  ),
  data.frame(
    structure = "arm", lambda_structure = 2L,
    lambda_log_ratio_sd_prior = prior_grid
  )
)
format_prior <- function(value) gsub("\\.", "p", sprintf("%.2f", value))
specifications$spec_id <- seq_len(nrow(specifications))
specifications$label <- ifelse(
  specifications$structure == "common",
  "common",
  paste0(
    specifications$structure,
    "-sd",
    vapply(
      specifications$lambda_log_ratio_sd_prior,
      format_prior,
      character(1)
    )
  )
)
specifications$seed <- 20260810L + specifications$spec_id * 100L
specifications$mode_status <- "pending"
specifications$sample_status <- "pending"
specifications$gq_status <- "pending"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(
  specifications,
  file.path(output_path, "lambda-specification-manifest.csv"),
  row.names = FALSE
)
print(specifications, row.names = FALSE)
