#!/usr/bin/env Rscript

# Summarize sample- and population-level social multipliers from the compact
# generated-quantities program. This script never fits or refits the model.
args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit)) substring(hit, nchar(name) + 2L) else default
}

gq_files <- strsplit(
  option_value("--gq-files", ""), ",", fixed = TRUE
)[[1L]]
distance_index <- as.integer(option_value("--distance-index", "1"))
distance_m <- as.numeric(option_value("--distance-m", "1250"))
output_path <- option_value(
  "--output-path", "temp-data/main-core-cluster-shock-population-1250"
)
table_path <- option_value(
  "--table-path",
  file.path(output_path, "main-core-cluster-shock-population-1250m.tex")
)

if (length(gq_files) != 4L || any(!nzchar(gq_files)) ||
    any(!file.exists(gq_files))) {
  stop("--gq-files must list four existing compact-GQ CSVs.", call. = FALSE)
}
if (!distance_index %in% 1:3) {
  stop("--distance-index must be 1, 2, or 3.", call. = FALSE)
}
if (!is.finite(distance_m) || distance_m < 0) {
  stop("--distance-m must be finite and nonnegative.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(readr)
  library(tidyr)
})

treatments <- c("Control", "Ink", "Calendar", "Bracelet")
object_catalog <- tibble::tribble(
  ~variable, ~estimand, ~threshold_probability,
  "core_compact_sm_sample_average", "Sample-average multiplier", TRUE,
  "core_compact_sm_population_average", "Population-average multiplier", TRUE,
  "core_compact_sm_sample_share_amplification",
    "Share of experimental communities amplifying", FALSE,
  "core_compact_sm_population_share_amplification",
    "Population probability of amplification", FALSE,
  "core_compact_sm_sample_share_mitigation",
    "Share of experimental communities mitigating", FALSE,
  "core_compact_sm_population_share_mitigation",
    "Population probability of mitigation", FALSE
)

variables <- unlist(lapply(object_catalog$variable, function(prefix) {
  sprintf("%s[%d,%d]", prefix, distance_index, seq_along(treatments))
}))
legacy_average <- sprintf(
  "core_compact_sm_rescaled[%d,%d]", distance_index, seq_along(treatments)
)

draws <- cmdstanr::read_cmdstan_csv(
  sort(gq_files), variables = c(variables, legacy_average)
)$generated_quantities |>
  posterior::as_draws_df() |>
  as.data.frame() |>
  mutate(draw = row_number())

# The new explicitly named sample average must reproduce the existing compact
# GQ estimand after applying the long-standing paper sign convention.
for (treatment_index in seq_along(treatments)) {
  sample_name <- sprintf(
    "core_compact_sm_sample_average[%d,%d]",
    distance_index, treatment_index
  )
  legacy_name <- sprintf(
    "core_compact_sm_rescaled[%d,%d]", distance_index, treatment_index
  )
  discrepancy <- max(abs(draws[[sample_name]] + draws[[legacy_name]]))
  if (!is.finite(discrepancy) || discrepancy > 1e-10) {
    stop(
      "New sample-average multiplier does not reproduce the existing GQ: ",
      sample_name, call. = FALSE
    )
  }
}

draw_level <- draws |>
  select(draw, all_of(variables)) |>
  pivot_longer(-draw, names_to = "stan_variable", values_to = "value") |>
  extract(
    stan_variable, c("variable", "row", "treatment_index"),
    "^(.+)\\[(\\d+),(\\d+)\\]$", convert = TRUE
  ) |>
  left_join(object_catalog, by = "variable") |>
  mutate(
    treatment = treatments[treatment_index],
    distance_m = distance_m
  ) |>
  select(draw, distance_m, estimand, treatment, value, threshold_probability)

summary <- draw_level |>
  group_by(distance_m, estimand, treatment, threshold_probability) |>
  summarise(
    median = median(value),
    conf_low = quantile(value, 0.025),
    conf_high = quantile(value, 0.975),
    posterior_probability_mitigation = if (first(threshold_probability)) {
      mean(value < 1)
    } else {
      NA_real_
    },
    draws = n(),
    .groups = "drop"
  )

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  draw_level |> select(-threshold_probability),
  file.path(output_path, "draws.csv")
)
readr::write_csv(summary, file.path(output_path, "summary.csv"))

format_probability <- function(x) {
  if (!is.finite(x)) return(NA_character_)
  if (x < 0.001) return("$<0.001$")
  if (x > 0.999) return("$>0.999$")
  sprintf("%.3f", x)
}
format_cell <- function(row) {
  if (row$threshold_probability) {
    sprintf(
      "\\makecell[c]{%.2f\\\\{[%.2f, %.2f]}\\\\%s}",
      row$median, row$conf_low, row$conf_high,
      format_probability(row$posterior_probability_mitigation)
    )
  } else {
    sprintf(
      "\\makecell[c]{%.2f\\\\{[%.2f, %.2f]}}",
      row$median, row$conf_low, row$conf_high
    )
  }
}

estimand_order <- object_catalog$estimand
treatment_order <- c("Control", "Calendar", "Ink", "Bracelet")
table_data <- summary |>
  mutate(
    estimand = factor(estimand, levels = estimand_order),
    treatment = factor(treatment, levels = treatment_order)
  ) |>
  arrange(estimand, treatment)

con <- file(table_path, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(c(
  "\\centering",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Estimand & Control & Calendar & Ink & Bracelet \\\\",
  "\\midrule"
), con)
for (estimand_name in estimand_order) {
  panel <- table_data |> filter(estimand == estimand_name)
  if (nrow(panel) != length(treatment_order)) next
  cells <- vapply(seq_len(nrow(panel)), function(i) {
    format_cell(panel[i, , drop = FALSE])
  }, character(1))
  writeLines(paste(c(estimand_name, cells), collapse = " & ") |> paste0(" \\\\"), con)
}
writeLines(c("\\bottomrule", "\\end{tabular}"), con)

message("Wrote cluster-shock multiplier summaries at ", distance_m, "m to ", output_path)
