#!/usr/bin/env Rscript

# Validate streamlined replacements for the four active robustness rows before
# any legacy posterior CSV is removed. This script reads only compact new fits
# and compact legacy generated quantities; it never loads the large legacy CSVs.

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit)) substring(hit, nchar(prefix) + 1L) else default
}

new_root <- value(
  "--new-root",
  "/project/akaring/takeup-data/data/stan_analysis_data/streamlined-active-robustness"
)
legacy_draws <- value(
  "--legacy-draws",
  "/project/akaring/takeup-data/candidate-multiplier-1250-20260825/multiplier-draws-1250.csv"
)
legacy_fit_root <- value(
  "--legacy-fit-root",
  "/project/akaring/takeup-data/data/stan_analysis_data"
)
corrected_legacy_root <- value(
  "--corrected-legacy-root",
  paste0(legacy_fit_root, "/streamlined-active-robustness-corrected-legacy")
)
output_path <- value("--output-path", file.path(new_root, "audit"))
median_tolerance <- as.numeric(value("--median-tolerance", "0.04"))
quantile_tolerance <- as.numeric(value("--quantile-tolerance", "0.07"))
probability_tolerance <- as.numeric(value("--probability-tolerance", "0.08"))

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(purrr)
  library(readr)
  library(tidyr)
})

catalog <- tibble::tribble(
  ~spec_id, ~schema, ~legacy_stem,
  "private-distance-community-image", "legacy",
    "dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS",
  "full-information", "legacy",
    "dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP",
  "exclude-dispersed", "core",
    "dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
  "second-order-observability", "core",
    "dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB"
)
treatments <- c("Control", "Ink", "Calendar", "Bracelet")

read_gq <- function(spec_id, schema, legacy_stem, root = new_root) {
  directory <- file.path(root, spec_id, "gq-1250")
  files <- list.files(directory, pattern = "[.]csv$", full.names = TRUE)
  files <- files[!grepl("diagnostics|manifest|audit", basename(files))]
  if (length(files) != 4L) {
    stop("Expected four compact GQ chains in ", directory, "; found ",
         length(files), call. = FALSE)
  }
  prefix <- if (schema == "core") {
    "core_compact_sm_rescaled"
  } else {
    "compact_sm_rescaled"
  }
  variables <- sprintf("%s[1,%d]", prefix, seq_along(treatments))
  draws <- cmdstanr::read_cmdstan_csv(sort(files), variables = variables)$generated_quantities |>
    posterior::as_draws_df() |>
    as.data.frame()
  draws |>
    transmute(draw = .data$.draw, across(all_of(variables))) |>
    pivot_longer(-.data$draw, names_to = "variable", values_to = "raw_multiplier") |>
    mutate(
      spec_id = spec_id,
      treatment = treatments[match(.data$variable, variables)],
      multiplier = -.data$raw_multiplier
    ) |>
    select(.data$spec_id, .data$treatment, .data$draw, .data$multiplier)
}

summarize_draws <- function(data, prefix) {
  data |>
    group_by(.data$spec_id, .data$treatment) |>
    summarise(
      "{prefix}_draws" := n(),
      "{prefix}_median" := median(.data$multiplier),
      "{prefix}_q025" := quantile(.data$multiplier, 0.025),
      "{prefix}_q975" := quantile(.data$multiplier, 0.975),
      "{prefix}_pr_lt_1" := mean(.data$multiplier < 1),
      .groups = "drop"
    )
}

if (!file.exists(legacy_draws)) stop("Missing legacy draw file: ", legacy_draws)
legacy <- readr::read_csv(legacy_draws, show_col_types = FALSE) |>
  filter(.data$table_id == "main", .data$spec_id %in% catalog$spec_id) |>
  select(.data$spec_id, .data$treatment, .data$draw, .data$multiplier)
new <- purrr::pmap_dfr(catalog, read_gq)

# The historical compact GQ hard-coded first-order belief coefficients in its
# multiplier block, including for the SOB model.  Compare the streamlined SOB
# result with the same old posterior evaluated by the corrected main-core GQ;
# otherwise this gate would test a known postprocessing bug rather than fit
# equivalence.
corrected_sob <- read_gq(
  "second-order-observability", "core",
  "dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB",
  root = corrected_legacy_root
)
legacy <- legacy |>
  filter(.data$spec_id != "second-order-observability") |>
  bind_rows(corrected_sob)
updated_draws <- readr::read_csv(legacy_draws, show_col_types = FALSE) |>
  filter(!(.data$table_id == "main" & .data$spec_id %in% catalog$spec_id)) |>
  bind_rows(new |>
    transmute(
      table_id = "main", .data$spec_id, .data$treatment, .data$draw,
      distance_m = 1250, .data$multiplier
    )) |>
  arrange(.data$table_id, .data$spec_id, .data$treatment, .data$draw)

comparison <- summarize_draws(legacy, "legacy") |>
  inner_join(summarize_draws(new, "new"), by = c("spec_id", "treatment")) |>
  mutate(
    median_difference = .data$new_median - .data$legacy_median,
    q025_difference = .data$new_q025 - .data$legacy_q025,
    q975_difference = .data$new_q975 - .data$legacy_q975,
    pr_lt_1_difference = .data$new_pr_lt_1 - .data$legacy_pr_lt_1,
    interval_overlap = pmax(
      0,
      pmin(.data$new_q975, .data$legacy_q975) -
        pmax(.data$new_q025, .data$legacy_q025)
    ) / pmax(
      .Machine$double.eps,
      pmax(.data$new_q975, .data$legacy_q975) -
        pmin(.data$new_q025, .data$legacy_q025)
    ),
    pass =
      abs(.data$median_difference) <= median_tolerance &
      abs(.data$q025_difference) <= quantile_tolerance &
      abs(.data$q975_difference) <= quantile_tolerance &
      abs(.data$pr_lt_1_difference) <= probability_tolerance &
      .data$interval_overlap >= 0.75
  ) |>
  mutate(
    reference_kind = if_else(
      .data$spec_id == "second-order-observability",
      "legacy_posterior_corrected_second_order_gq",
      "historical_compact_gq"
    ),
    .after = .data$treatment
  )

diagnose_spec <- function(spec_id, schema, legacy_stem) {
  fit_dir <- file.path(new_root, spec_id, "fits")
  files <- list.files(
    fit_dir, pattern = "-slim-chain[1-4]-1[.]csv$", full.names = TRUE
  )
  if (length(files) != 4L) stop("Expected four fit chains for ", spec_id)
  fit <- cmdstanr::read_cmdstan_csv(sort(files))$post_warmup_draws
  fit_summary <- posterior::summarise_draws(fit)
  diagnostics <- cmdstanr::read_cmdstan_csv(
    sort(files)
  )$post_warmup_sampler_diagnostics
  diagnostic_df <- posterior::as_draws_df(diagnostics)
  energy <- posterior::as_draws_array(diagnostics[, , "energy__"])
  ebfmi <- vapply(seq_len(dim(energy)[2]), function(chain) {
    values <- energy[, chain, 1]
    mean(diff(values)^2) / stats::var(values)
  }, numeric(1))
  tibble(
    spec_id = spec_id,
    replacement_method = if (schema == "core") "main_core_refit" else
      "lossless_column_subset",
    chains = length(files),
    retained_draws = posterior::ndraws(fit),
    max_rhat = max(fit_summary$rhat, na.rm = TRUE),
    min_bulk_ess = min(fit_summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(fit_summary$ess_tail, na.rm = TRUE),
    divergences = sum(diagnostic_df$divergent__),
    max_treedepth_hits = sum(diagnostic_df$treedepth__ >= 12),
    min_ebfmi = min(ebfmi),
    fit_bytes = sum(file.info(files)$size)
  )
}
diagnostics <- purrr::pmap_dfr(catalog, diagnose_spec)

deletion_manifest <- purrr::pmap_dfr(catalog, function(spec_id, schema, legacy_stem) {
  paths <- file.path(legacy_fit_root, paste0(legacy_stem, "-", 1:4, ".csv"))
  info <- file.info(paths)
  tibble(
    spec_id = spec_id,
    path = paths,
    exists = file.exists(paths),
    bytes = info$size,
    replacement_root = file.path(new_root, spec_id),
    delete_approved_by_audit = all(comparison$pass[comparison$spec_id == spec_id]) &&
      all(file.exists(paths)) &&
      (schema != "core" || (
        diagnostics$divergences[diagnostics$spec_id == spec_id] == 0 &&
        diagnostics$max_treedepth_hits[diagnostics$spec_id == spec_id] == 0
      ))
  )
})

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(new, file.path(output_path, "streamlined-multiplier-draws-1250.csv"))
readr::write_csv(updated_draws, file.path(output_path, "updated-multiplier-draws-1250.csv"))
readr::write_csv(comparison, file.path(output_path, "legacy-streamlined-comparison.csv"))
readr::write_csv(diagnostics, file.path(output_path, "streamlined-fit-diagnostics.csv"))
readr::write_csv(deletion_manifest, file.path(output_path, "legacy-deletion-manifest.csv"))

if (any(!comparison$pass)) {
  stop("At least one result-equivalence check failed; legacy deletion is forbidden.")
}
core_diagnostics <- diagnostics |>
  filter(.data$replacement_method == "main_core_refit")
if (any(core_diagnostics$divergences > 0 |
        core_diagnostics$max_treedepth_hits > 0)) {
  stop("A new main-core refit has hard HMC failures; legacy deletion is forbidden.")
}
if (any(!deletion_manifest$delete_approved_by_audit)) {
  stop("Deletion manifest is incomplete or unapproved.")
}
message("All streamlined equivalence and hard-diagnostic gates passed.")
