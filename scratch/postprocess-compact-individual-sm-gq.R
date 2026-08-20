#!/usr/bin/env Rscript

# Convert compact individual-model GQ CSVs to the same RDS schema produced by
# quick_roc_postprocess.R --sm.

args <- commandArgs(trailingOnly = TRUE)

option_value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

model <- option_value("--model")
input_path <- option_value(
  "--input-path",
  "/project/akaring/takeup-data/temp-data/compact-sm-gq"
)
output_path <- option_value(
  "--output-path",
  "/project/akaring/takeup-data/temp-data"
)

allowed_models <- c(
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS",
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP"
)
if (is.null(model) || !model %in% allowed_models) {
  stop(
    "--model must be one of: ",
    paste(allowed_models, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(tidyr)
  library(tidybayes)
})

csv_files <- file.path(
  input_path,
  sprintf("compact_sm_fit105_%s-%d.csv", model, 1:4)
)
if (!all(file.exists(csv_files))) {
  stop(
    "Missing compact GQ chain(s): ",
    paste(basename(csv_files[!file.exists(csv_files)]), collapse = ", "),
    call. = FALSE
  )
}

fit <- cmdstanr::read_cmdstan_csv(csv_files)$generated_quantities
draws <- spread_rvars(
  fit,
  compact_sm[compact_dist_idx, treat_idx],
  compact_dist_beta_v[compact_dist_idx, treat_idx],
  compact_sm_delta_part[compact_dist_idx, treat_idx],
  compact_sm_mu_part[compact_dist_idx, treat_idx],
  compact_sm_rescaled[compact_dist_idx, treat_idx],
  compact_sm_delta_part_rescaled[compact_dist_idx, treat_idx],
  compact_sm_mu_part_rescaled[compact_dist_idx, treat_idx]
) |>
  mutate(
    model = model,
    fit_version = 105L,
    fit_type = "fit",
    roc_distance = c(500, 1500, 2500)[compact_dist_idx],
    treatment = factor(
      c("Control", "Ink", "Calendar", "Bracelet")[treat_idx],
      levels = c("Control", "Ink", "Calendar", "Bracelet")
    )
  ) |>
  select(
    model,
    fit_version,
    fit_type,
    treatment,
    roc_distance,
    sm = compact_sm,
    dist_beta_v = compact_dist_beta_v,
    sm_delta_part = compact_sm_delta_part,
    sm_mu_part = compact_sm_mu_part,
    sm_rescaled = compact_sm_rescaled,
    sm_delta_part_rescaled = compact_sm_delta_part_rescaled,
    sm_mu_part_rescaled = compact_sm_mu_part_rescaled
  ) |>
  pivot_longer(where(is_rvar), names_to = "variable")

expected_rows <- 4L * 3L * 7L
if (nrow(draws) != expected_rows) {
  stop("Compact draw grid has ", nrow(draws), " rows; expected ", expected_rows, ".")
}

summary <- draws |>
  median_qi(value) |>
  tidybayes::to_broom_names()

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
draw_path <- file.path(
  output_path,
  sprintf(
    "rvar_processed_dist_fit105_sm_draws_%s_1-4.rds",
    model
  )
)
summary_path <- file.path(
  output_path,
  sprintf(
    "rvar_processed_dist_fit105_sm_summ_%s_1-4.rds",
    model
  )
)
saveRDS(draws, draw_path)
saveRDS(summary, summary_path)

message("Wrote compact draw object: ", draw_path)
message("Wrote compact summary object: ", summary_path)
