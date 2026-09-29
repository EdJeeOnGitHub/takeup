#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
})
source("R/reduced-form/functions.R")
# Formatting is irrelevant to the significance thresholds under test.
linebreak <- function(x, align) x
p <- c(0.0096, 0.0496, 0.0996, 0.1004)
fixture <- tibble(
  assigned_treatment = "bracelet",
  assigned_dist_group = c("combined", "close", "far", "far - close"),
  estimate = qnorm(1 - p / 2), std_error = 1,
  conf.low = estimate - 1.96, conf.high = estimate + 1.96,
  pval = round_pval(p), oneside_pval = round_pval(p / 2),
  show_pval_only = FALSE, n_obs_line = FALSE
)
result <- suppressWarnings(prep_tbl(fixture, stat = "std.error", stars = TRUE))
values <- unlist(result[1, c("Combined", "Close", "Far", "Far - Close")])
star_counts <- str_count(values, fixed("*"))
stopifnot(identical(unname(star_counts), c(3L, 2L, 1L, 0L)))
cat("Table stars match the manuscript thresholds using unrounded p-values.\n")
