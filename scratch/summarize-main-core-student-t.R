#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

baseline_gq_dir <- option_value("--baseline-gq-dir")
student_fit_dir <- option_value("--student-fit-dir")
student_gq_dir <- option_value("--student-gq-dir")
output_path <- option_value(
  "--output-path", "temp-data/main-core-student-t-robustness"
)
if (any(vapply(
  c(baseline_gq_dir, student_fit_dir, student_gq_dir),
  function(path) is.null(path) || !dir.exists(path), logical(1)
))) stop("Existing baseline GQ, Student-t fit, and Student-t GQ directories are required.",
        call. = FALSE)

csvs <- function(path, pattern = "\\.csv$") {
  sort(list.files(path, pattern = pattern, full.names = TRUE))
}
baseline_gq_csvs <- csvs(baseline_gq_dir)
student_fit_csvs <- csvs(student_fit_dir, "student-t5-chain[1-4]-.*\\.csv$")
student_gq_csvs <- csvs(student_gq_dir)
if (length(baseline_gq_csvs) != 4L || length(student_fit_csvs) != 4L ||
    length(student_gq_csvs) != 4L) {
  stop("Expected four CSVs for each fit/GQ input.", call. = FALSE)
}

baseline <- as_draws_df(read_cmdstan_csv(baseline_gq_csvs)$generated_quantities)
student_gq <- as_draws_df(read_cmdstan_csv(student_gq_csvs)$generated_quantities)
student_fit <- read_cmdstan_csv(student_fit_csvs)
student_summary <- summarise_draws(student_fit$post_warmup_draws)
student_sampler <- as_draws_df(student_fit$post_warmup_sampler_diagnostics)

treatments <- c("Control", "Ink", "Calendar", "Bracelet")
distances <- c("500m", "1500m", "2500m")
summarize_value <- function(value) c(
  median = median(value),
  q025 = unname(quantile(value, 0.025)),
  q975 = unname(quantile(value, 0.975))
)
rows <- list()
for (specification in c("Gaussian", "Student-t(5)")) {
  draws <- if (specification == "Gaussian") baseline else student_gq
  for (treatment_index in 1:4) {
    for (distance_index in 1:3) {
      value <- -draws[[sprintf(
        "core_compact_sm_rescaled[%d,%d]", distance_index, treatment_index
      )]]
      summary <- summarize_value(value)
      rows[[length(rows) + 1L]] <- data.frame(
        specification = specification,
        treatment = treatments[treatment_index],
        distance = distances[distance_index],
        median = summary["median"], q025 = summary["q025"],
        q975 = summary["q975"]
      )
    }
  }
}
results <- bind_rows(rows)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output_path, "student-t-multiplier-summary.csv"),
          row.names = FALSE)
diagnostics <- data.frame(
  max_rhat = max(student_summary$rhat, na.rm = TRUE),
  min_ess_bulk = min(student_summary$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(student_summary$ess_tail, na.rm = TRUE),
  divergences = sum(student_sampler$divergent__, na.rm = TRUE),
  max_treedepth = sum(student_sampler$treedepth__ >= 12, na.rm = TRUE)
)
write.csv(diagnostics, file.path(output_path, "student-t-fit-diagnostics.csv"),
          row.names = FALSE)

fmt <- function(x) sprintf("%.2f", x)
ci <- function(lo, hi) paste0("(", fmt(lo), ", ", fmt(hi), ")")
lines <- c(
  "\\begin{tabular}{lccc}", "\\toprule",
  "Treatment & 500m & 1500m & 2500m \\\\",
  "\\midrule"
)
for (specification in c("Gaussian", "Student-t(5)")) {
  block <- results |> filter(.data$specification == specification)
  panel <- if (specification == "Gaussian") {
    "Panel A: Gaussian intrinsic motivation"
  } else {
    "Panel B: Student-$t(5)$ intrinsic motivation"
  }
  lines <- c(lines, paste0("\\multicolumn{4}{l}{\\textit{", panel, "}} \\\\"))
  for (treatment_name in c("Bracelet", "Calendar", "Ink", "Control")) {
    cells <- block |>
      filter(.data$treatment == .env$treatment_name) |>
      arrange(match(.data$distance, distances))
    lines <- c(
      lines,
      paste0(
        treatment_name, " & ",
        paste(fmt(cells$median), collapse = " & "), " \\\\"
      ),
      paste0(" & ", paste(ci(cells$q025, cells$q975), collapse = " & "),
             " \\\\ ")
    )
  }
  if (specification == "Gaussian") lines <- c(lines, "\\addlinespace")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(output_path, "main-core-student-t-multipliers.tex"))
print(diagnostics)
