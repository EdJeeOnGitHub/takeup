#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
option_value <- function(name, default = NULL) main_core_option_value(args, name, default)
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

specifications <- data.frame(
  id = c("baseline", "conditional", "unconditional"),
  label = c(
    "Perfect classification (baseline)",
    "Asymmetric reports, conditional on recognition",
    "Asymmetric reports, unrecognized as null signal"
  ),
  fit_dir = c(
    option_value("--baseline-fit-dir"),
    option_value("--conditional-fit-dir"),
    option_value("--unconditional-fit-dir")
  ),
  gq_dir = c(
    option_value("--baseline-gq-dir"),
    option_value("--conditional-gq-dir"),
    option_value("--unconditional-gq-dir")
  ),
  stringsAsFactors = FALSE
)
output_path <- option_value(
  "--output-path", "temp-data/main-core-asymmetric-observability"
)
if (any(is.na(specifications$fit_dir)) || any(is.na(specifications$gq_dir)) ||
    any(!dir.exists(c(specifications$fit_dir, specifications$gq_dir)))) {
  stop("All three fit and compact-GQ directories are required.", call. = FALSE)
}

csvs <- function(path) sort(list.files(path, "\\.csv$", full.names = TRUE))
summarize_value <- function(value) c(
  median = median(value, na.rm = TRUE),
  lower = unname(quantile(value, 0.025, na.rm = TRUE)),
  upper = unname(quantile(value, 0.975, na.rm = TRUE))
)
treatments <- c("Control", "Ink", "Calendar", "Bracelet")
distance_labels <- c("500m", "1500m", "2500m")
level_labels <- c("Combined", "Close", "Far")
multiplier_rows <- list()
ate_rows <- list()
response_rows <- list()
diagnostic_rows <- list()

for (specification_index in seq_len(nrow(specifications))) {
  specification <- specifications[specification_index, ]
  fit_files <- csvs(specification$fit_dir)
  gq_files <- csvs(specification$gq_dir)
  if (length(fit_files) != 4L || length(gq_files) != 4L) {
    stop("Expected four fit and GQ chains for ", specification$id, call. = FALSE)
  }
  fit <- read_cmdstan_csv(fit_files)
  gq <- as_draws_df(read_cmdstan_csv(gq_files)$generated_quantities)
  fit_summary <- summarise_draws(fit$post_warmup_draws)
  sampler <- as_draws_df(fit$post_warmup_sampler_diagnostics)
  diagnostic_rows[[specification_index]] <- data.frame(
    specification = specification$label,
    max_rhat = max(fit_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(fit_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(fit_summary$ess_tail, na.rm = TRUE),
    divergences = sum(sampler$divergent__, na.rm = TRUE),
    max_treedepth = sum(sampler$treedepth__ >= 12, na.rm = TRUE)
  )
  for (treatment_index in seq_along(treatments)) {
    for (distance_index in seq_along(distance_labels)) {
      multiplier <- -gq[[sprintf(
        "core_compact_sm_rescaled[%d,%d]", distance_index, treatment_index
      )]]
      multiplier_summary <- summarize_value(multiplier)
      multiplier_rows[[length(multiplier_rows) + 1L]] <- data.frame(
        specification = specification$label,
        treatment = treatments[treatment_index],
        distance = distance_labels[distance_index],
        t(multiplier_summary), check.names = FALSE
      )
      if (specification$id != "baseline") {
        response <- lapply(c(
          "recognition_nontaker", "recognition_taker", "sensitivity",
          "specificity", "information_factor"
        ), function(quantity) summarize_value(gq[[sprintf(
          "core_compact_%s[%d,%d]", quantity, distance_index, treatment_index
        )]]))
        names(response) <- c(
          "recognition_nontaker", "recognition_taker", "sensitivity",
          "specificity", "information_factor"
        )
        for (quantity in names(response)) {
          response_rows[[length(response_rows) + 1L]] <- data.frame(
            specification = specification$label,
            treatment = treatments[treatment_index],
            distance = distance_labels[distance_index], quantity = quantity,
            t(response[[quantity]]), check.names = FALSE
          )
        }
      }
    }
  }
  for (treatment_index in 2:4) {
    level_draws <- vapply(seq_along(level_labels), function(level_index) {
      gq[[sprintf("core_compact_takeup_level[%d,%d]", level_index, treatment_index)]] -
        gq[[sprintf("core_compact_takeup_level[%d,1]", level_index)]]
    }, numeric(nrow(gq)))
    colnames(level_draws) <- level_labels
    effect_draws <- cbind(
      level_draws,
      `Far - Close` = level_draws[, "Far"] - level_draws[, "Close"]
    )
    for (effect in colnames(effect_draws)) {
      summary <- summarize_value(effect_draws[, effect])
      ate_rows[[length(ate_rows) + 1L]] <- data.frame(
        specification = specification$label,
        treatment = treatments[treatment_index], effect = effect,
        t(summary), check.names = FALSE
      )
    }
  }
}

multiplier <- do.call(rbind, multiplier_rows)
ate <- do.call(rbind, ate_rows)
response <- do.call(rbind, response_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(multiplier, file.path(output_path, "asymmetric-multiplier-summary.csv"), row.names = FALSE)
write.csv(ate, file.path(output_path, "asymmetric-ate-summary.csv"), row.names = FALSE)
write.csv(response, file.path(output_path, "asymmetric-response-summary.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_path, "asymmetric-fit-diagnostics.csv"), row.names = FALSE)

fmt <- function(x) sprintf("%.2f", x)
interval <- function(lo, hi) paste0("(", fmt(lo), ", ", fmt(hi), ")")
panel_table <- function(data, columns, row_names, file, caption_label) {
  lines <- c(
    paste0("\\begin{tabular}{l", paste(rep("c", length(columns)), collapse = ""), "}"),
    "\\toprule",
    paste(c("Treatment", columns), collapse = " & ") |> paste0(" \\\\"),
    "\\midrule"
  )
  for (specification_index in seq_len(nrow(specifications))) {
    block <- data[data$specification == specifications$label[specification_index], ]
    if (specification_index > 1L) lines <- c(lines, "\\addlinespace")
    lines <- c(lines, paste0(
      "\\multicolumn{", length(columns) + 1L, "}{l}{\\textit{Panel ",
      LETTERS[specification_index], ": ", specifications$label[specification_index],
      "}} \\\\"
    ))
    for (row_name in row_names) {
      cells <- block[block$treatment == row_name, ]
      order_index <- match(columns, if ("distance" %in% names(cells)) cells$distance else cells$effect)
      cells <- cells[order_index, ]
      lines <- c(
        lines,
        paste0(row_name, " & ", paste(fmt(cells$median), collapse = " & "), " \\\\"),
        paste0(" & ", paste(interval(cells$lower, cells$upper), collapse = " & "), " \\\\ ")
      )
    }
  }
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}",
               paste0("% ", caption_label)), file.path(output_path, file))
}
panel_table(
  multiplier, distance_labels, rev(treatments),
  "main-core-asymmetric-multipliers.tex", "Asymmetric observability multipliers"
)
panel_table(
  ate, c(level_labels, "Far - Close"), c("Bracelet", "Calendar", "Ink"),
  "main-core-asymmetric-ates.tex", "Asymmetric observability treatment effects"
)
print(diagnostics)
