#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}
suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(purrr)
  library(rlang)
})

root <- option_value(
  "--root",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder"
)
f0_fit <- option_value(
  "--f0-fit-dir",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-asym-conditional-production"
)
f0_gq <- option_value("--f0-gq-dir", file.path(f0_fit, "compact-gq"))
output_path <- option_value("--output-path", file.path(root, "summary"))
workspace <- option_value(
  "--workspace",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-asym-input/dist_fit105.RData"
)
summary_data <- prepare_main_core_data(
  workspace, "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
  observation_model = 1L
)
finite_distance_change <- summary_data$roc_distances[26L] -
  summary_data$roc_distances[6L]

specifications <- data.frame(
  id = c("f0", "f1", "f2", "f3", "u3"),
  label = c(
    "F0: Fully flexible conditional channel",
    "F1: Restricted multinomial conditional channel",
    "F2: Recognition-conditioned multinomial channel",
    "F3: Two-stage conditional channel",
    "U3: Two-stage unconditional channel"
  ),
  observation_model = c(1L, 1L, 1L, 1L, 2L),
  recognition_structure = c(0L, 0L, 2L, 2L, 1L),
  report_structure = c(0L, 1L, 1L, 2L, 2L),
  measurement_parameters = c(48L, 36L, 20L, 19L, 21L),
  fit_dir = c(f0_fit, file.path(root, c("f1", "f2", "f3", "u3"))),
  gq_dir = c(f0_gq, file.path(root, c("f1", "f2", "f3", "u3"), "gq")),
  stringsAsFactors = FALSE
)

csvs <- function(path, fit = FALSE) {
  files <- sort(list.files(path, "[.]csv$", full.names = TRUE))
  files <- files[!grepl("peer-link-audit|profile|status", basename(files))]
  if (fit) files <- files[grepl("chain[1-4]-1[.]csv$", basename(files))]
  files[file.info(files)$size > 10000]
}
summarize_value <- function(value) c(
  median = median(value, na.rm = TRUE),
  lower = unname(quantile(value, 0.025, na.rm = TRUE)),
  upper = unname(quantile(value, 0.975, na.rm = TRUE)),
  finite_share = mean(is.finite(value))
)
treatments <- c("Control", "Ink", "Calendar", "Bracelet")
distances <- c("500m", "1500m", "2500m")
levels <- c("Combined", "Close", "Far")
multiplier_rows <- ate_rows <- response_rows <- diagnostic_rows <- list()

for (specification_index in seq_len(nrow(specifications))) {
  specification <- specifications[specification_index, ]
  fit_files <- csvs(specification$fit_dir, fit = TRUE)
  gq_files <- csvs(specification$gq_dir)
  if (length(fit_files) != 4L || length(gq_files) != 4L) {
    stop("Expected four fit and GQ files for ", specification$id, call. = FALSE)
  }
  fit <- read_cmdstan_csv(fit_files)
  fit_summary <- summarise_draws(fit$post_warmup_draws)
  sampler <- as_draws_df(fit$post_warmup_sampler_diagnostics)
  diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
    id = specification$id, specification = specification$label,
    max_rhat = max(fit_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(fit_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(fit_summary$ess_tail, na.rm = TRUE),
    divergences = sum(sampler$divergent__, na.rm = TRUE),
    max_treedepth = sum(sampler$treedepth__ >= 12, na.rm = TRUE)
  )
  gq <- as_draws_df(read_cmdstan_csv(gq_files)$generated_quantities)
  for (treatment_index in seq_along(treatments)) {
    for (distance_index in seq_along(distances)) {
      value <- -gq[[sprintf(
        "core_compact_sm_rescaled[%d,%d]", distance_index, treatment_index
      )]]
      multiplier_rows[[length(multiplier_rows) + 1L]] <- data.frame(
        id = specification$id, specification = specification$label,
        treatment = treatments[treatment_index], estimand = "Point",
        distance = distances[distance_index],
        t(summarize_value(value)), check.names = FALSE
      )
    }
    finite_name <- sprintf("core_compact_sm_finite[%d]", treatment_index)
    finite <- if (finite_name %in% names(gq)) {
      gq[[finite_name]]
    } else {
      numerator <- gq[[sprintf("core_compact_cutoff[3,%d]", treatment_index)]] -
        gq[[sprintf("core_compact_cutoff[1,%d]", treatment_index)]]
      dist_beta <- as_draws_df(fit$post_warmup_draws)[["dist_beta_v[1]"]]
      numerator / (dist_beta * finite_distance_change)
    }
    multiplier_rows[[length(multiplier_rows) + 1L]] <- data.frame(
      id = specification$id, specification = specification$label,
      treatment = treatments[treatment_index], estimand = "Finite",
      distance = "500--2500m", t(summarize_value(finite)), check.names = FALSE
    )
  }
  for (treatment_index in 2:4) {
    level_draws <- vapply(seq_along(levels), function(level_index) {
      gq[[sprintf(
        "core_compact_takeup_level[%d,%d]", level_index, treatment_index
      )]] - gq[[sprintf("core_compact_takeup_level[%d,1]", level_index)]]
    }, numeric(nrow(gq)))
    colnames(level_draws) <- levels
    effect_draws <- cbind(
      level_draws,
      `Far - Close` = level_draws[, "Far"] - level_draws[, "Close"]
    )
    for (effect in colnames(effect_draws)) {
      ate_rows[[length(ate_rows) + 1L]] <- data.frame(
        id = specification$id, specification = specification$label,
        treatment = treatments[treatment_index], effect = effect,
        t(summarize_value(effect_draws[, effect])), check.names = FALSE
      )
    }
  }
  if (specification$id != "f0") {
    for (quantity in c(
      "recognition_nontaker", "recognition_taker", "sensitivity",
      "specificity", "information_factor"
    )) for (treatment_index in seq_along(treatments)) {
      for (distance_index in seq_along(distances)) {
        value <- gq[[sprintf(
          "core_compact_%s[%d,%d]", quantity, distance_index, treatment_index
        )]]
        response_rows[[length(response_rows) + 1L]] <- data.frame(
          id = specification$id, specification = specification$label,
          treatment = treatments[treatment_index],
          distance = distances[distance_index], quantity = quantity,
          t(summarize_value(value)), check.names = FALSE
        )
      }
    }
  }
}

multiplier <- do.call(rbind, multiplier_rows)
ate <- do.call(rbind, ate_rows)
response <- do.call(rbind, response_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(specifications, file.path(output_path, "ladder-manifest.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_path, "ladder-diagnostics.csv"), row.names = FALSE)
write.csv(multiplier, file.path(output_path, "ladder-multipliers.csv"), row.names = FALSE)
write.csv(ate, file.path(output_path, "ladder-ates.csv"), row.names = FALSE)
write.csv(response, file.path(output_path, "ladder-response-channel.csv"), row.names = FALSE)

fmt <- function(x) sprintf("%.2f", x)
ci <- function(lower, upper) paste0("(", fmt(lower), ", ", fmt(upper), ")")
columns <- c("500m", "1500m", "2500m", "500--2500m")
lines <- c(
  "\\begin{tabular}{lcccc}", "\\toprule",
  "Treatment & Point: 500m & Point: 1500m & Point: 2500m & Finite: 500--2500m \\\\",
  "\\midrule"
)
for (specification_index in seq_len(nrow(specifications))) {
  block <- multiplier[multiplier$id == specifications$id[specification_index], ]
  if (specification_index > 1L) lines <- c(lines, "\\addlinespace")
  lines <- c(lines, paste0(
    "\\multicolumn{5}{l}{\\textit{Panel ", LETTERS[specification_index],
    ": ", specifications$label[specification_index], "}} \\\\"
  ))
  for (treatment in treatments) {
    cells <- block[block$treatment == treatment, ]
    cells <- cells[match(columns, cells$distance), ]
    lines <- c(
      lines,
      paste0(treatment, " & ", paste(fmt(cells$median), collapse = " & "), " \\\\"),
      paste0(" & ", paste(ci(cells$lower, cells$upper), collapse = " & "), " \\\\ ")
    )
  }
}
writeLines(
  c(lines, "\\bottomrule", "\\end{tabular}"),
  file.path(output_path, "main-core-observability-ladder-multipliers.tex")
)

# A compact forest plot makes the local-versus-finite comparison visible
# without requiring the paper's plotting stack.
pdf(
  file.path(output_path, "main-core-observability-ladder-multipliers.pdf"),
  width = 8, height = 8
)
par(mfrow = c(2, 2), mar = c(4, 5, 3, 1), las = 1)
for (treatment in treatments) {
  block <- multiplier[
    multiplier$treatment == treatment &
      multiplier$distance %in% c("1500m", "500--2500m"),
  ]
  point <- block[block$distance == "1500m", ][
    match(specifications$id, block$id[block$distance == "1500m"]),
  ]
  finite <- block[block$distance == "500--2500m", ][
    match(specifications$id, block$id[block$distance == "500--2500m"]),
  ]
  y <- rev(seq_len(nrow(specifications)))
  x_range <- range(c(point$lower, point$upper, finite$lower, finite$upper, 1))
  plot(
    NA, xlim = x_range, ylim = c(0.5, nrow(specifications) + 0.5),
    yaxt = "n", ylab = "", xlab = "Social multiplier", main = treatment
  )
  axis(2, at = y, labels = toupper(specifications$id))
  abline(v = 1, col = "grey65", lty = 2)
  segments(point$lower, y + 0.10, point$upper, y + 0.10, lwd = 1.5)
  points(point$median, y + 0.10, pch = 16)
  segments(finite$lower, y - 0.10, finite$upper, y - 0.10, lwd = 1.5)
  points(finite$median, y - 0.10, pch = 1)
  if (treatment == treatments[1]) {
    legend(
      "topright", c("Point at 1500m", "Finite 500--2500m"),
      pch = c(16, 1), bty = "n", cex = 0.8
    )
  }
}
dev.off()
print(diagnostics)
