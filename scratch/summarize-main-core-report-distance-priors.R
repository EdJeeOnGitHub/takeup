#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
source("scratch/main-core-multiplier-contrasts.R")
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
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors"
)
f0_fit <- option_value(
  "--f0-fit-dir",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-asym-conditional-production"
)
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

grouped_fit <- file.path(root, "design-pooled")
grouped_gq <- file.path(
  grouped_fit,
  if (dir.exists(file.path(grouped_fit, "gq-bracketed"))) {
    "gq-bracketed"
  } else {
    "gq"
  }
)
specifications <- data.frame(
  id = c("f0", "hierarchical", "tight", "design-pooled"),
  label = c(
    "Unrestricted multinomial, N(0, 0.25)",
    "Hierarchical multinomial, half-N(0, 0.25) scales",
    "Tighter multinomial, N(0, 0.10)",
    "Design-pooled multinomial"
  ),
  hierarchical = c(0L, 1L, 0L, 2L),
  prior_scale = c(0.25, 0.25, 0.10, 0.25),
  fit_dir = c(
    f0_fit, file.path(root, c("hierarchical", "tight")), grouped_fit
  ),
  gq_dir = c(
    file.path(f0_fit, "compact-gq"),
    file.path(root, c("hierarchical", "tight"), "gq"),
    grouped_gq
  ),
  stringsAsFactors = FALSE
)
specifications <- specifications[
  dir.exists(specifications$fit_dir) & dir.exists(specifications$gq_dir),
]

csvs <- function(path, fit = FALSE) {
  value <- sort(list.files(path, "[.]csv$", full.names = TRUE))
  value <- value[!grepl("peer-link-audit|profile|status", basename(value))]
  if (fit) value <- value[grepl("chain[1-4]-1[.]csv$", basename(value))]
  value[file.info(value)$size > 10000]
}
summarize_value <- function(value) {
  finite <- value[is.finite(value)]
  c(
    median = median(finite),
    lower = unname(quantile(finite, 0.025)),
    upper = unname(quantile(finite, 0.975)),
    probability_negative = mean(finite < 0),
    probability_below_one = mean(finite < 1),
    finite_share = length(finite) / length(value)
  )
}
treatments <- c("Control", "Ink", "Calendar", "Bracelet")
distances <- c("500m", "1500m", "2500m")
levels <- c("Combined", "Close", "Far")
multiplier_rows <- ate_rows <- diagnostic_rows <- scale_rows <- list()
contrast_rows <- list()

for (specification_index in seq_len(nrow(specifications))) {
  specification <- specifications[specification_index, ]
  fit_files <- csvs(specification$fit_dir, fit = TRUE)
  gq_files <- csvs(specification$gq_dir)
  if (length(fit_files) != 4L || length(gq_files) != 4L) {
    stop("Expected four fit and GQ files for ", specification$id, call. = FALSE)
  }
  fit <- read_cmdstan_csv(fit_files)
  fit_draws <- as_draws_df(fit$post_warmup_draws)
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
  if (specification$hierarchical == 1L) {
    for (truth in 1:2) for (category in 1:2) {
      parameter <- sprintf("core_report_arm_dist_sd[%d,%d]", truth, category)
      scale_rows[[length(scale_rows) + 1L]] <- data.frame(
        id = specification$id, truth = c("Nontaker", "Taker")[truth],
        report = c("Yes", "No")[category],
        t(summarize_value(fit_draws[[parameter]])), check.names = FALSE
      )
    }
  } else if (specification$hierarchical == 2L) {
    scale_rows[[length(scale_rows) + 1L]] <- data.frame(
      id = specification$id, truth = "All", report = "Within-group SD",
      t(summarize_value(fit_draws[["core_report_within_dist_sd[1]"]])),
      check.names = FALSE
    )
  }
  gq <- as_draws_df(read_cmdstan_csv(gq_files)$generated_quantities)
  point_multiplier <- matrix(
    NA_real_, nrow = nrow(gq), ncol = length(treatments) * length(distances)
  )
  dim(point_multiplier) <- c(nrow(gq), length(distances), length(treatments))
  finite_multiplier <- matrix(
    NA_real_, nrow = nrow(gq), ncol = length(treatments)
  )
  for (treatment_index in seq_along(treatments)) {
    for (distance_index in seq_along(distances)) {
      value <- -gq[[sprintf(
        "core_compact_sm_rescaled[%d,%d]", distance_index, treatment_index
      )]]
      point_multiplier[, distance_index, treatment_index] <- value
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
      numerator / (fit_draws[["dist_beta_v[1]"]] * finite_distance_change)
    }
    finite_multiplier[, treatment_index] <- finite
    multiplier_rows[[length(multiplier_rows) + 1L]] <- data.frame(
      id = specification$id, specification = specification$label,
      treatment = treatments[treatment_index], estimand = "Finite",
      distance = "500--2500m", t(summarize_value(finite)), check.names = FALSE
    )
  }
  primary_by_distance <- matrix(NA_real_, nrow(gq), length(distances))
  for (distance_index in seq_along(distances)) {
    contrasts <- main_core_multiplier_contrasts(
      point_multiplier[, distance_index, ]
    )
    primary_by_distance[, distance_index] <-
      contrasts[, "No Signal - Any Signal"]
    for (contrast_name in colnames(contrasts)) {
      contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
        id = specification$id, specification = specification$label,
        contrast = contrast_name, estimand = "Point",
        distance = distances[distance_index],
        t(main_core_summarize_contrast(contrasts[, contrast_name])),
        check.names = FALSE
      )
    }
  }
  finite_contrasts <- main_core_multiplier_contrasts(finite_multiplier)
  for (contrast_name in colnames(finite_contrasts)) {
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      id = specification$id, specification = specification$label,
      contrast = contrast_name, estimand = "Finite",
      distance = "500--2500m",
      t(main_core_summarize_contrast(finite_contrasts[, contrast_name])),
      check.names = FALSE
    )
  }
  minimum_primary <- apply(primary_by_distance, 1L, min, na.rm = FALSE)
  contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
    id = specification$id, specification = specification$label,
    contrast = "No Signal - Any Signal", estimand = "Grid minimum",
    distance = "Minimum", t(main_core_summarize_contrast(minimum_primary)),
    check.names = FALSE
  )
  for (treatment_index in 2:4) {
    level_draws <- vapply(seq_along(levels), function(level_index) {
      gq[[sprintf(
        "core_compact_takeup_level[%d,%d]", level_index, treatment_index
      )]] - gq[[sprintf("core_compact_takeup_level[%d,1]", level_index)]]
    }, numeric(nrow(gq)))
    colnames(level_draws) <- levels
    effects <- cbind(
      level_draws,
      `Far - Close` = level_draws[, "Far"] - level_draws[, "Close"]
    )
    for (effect in colnames(effects)) {
      ate_rows[[length(ate_rows) + 1L]] <- data.frame(
        id = specification$id, specification = specification$label,
        treatment = treatments[treatment_index], effect = effect,
        t(summarize_value(effects[, effect])), check.names = FALSE
      )
    }
  }
}

multiplier <- do.call(rbind, multiplier_rows)
ate <- do.call(rbind, ate_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
scales <- if (length(scale_rows)) do.call(rbind, scale_rows) else data.frame()
contrasts <- do.call(rbind, contrast_rows)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(specifications, file.path(output_path, "prior-specification-manifest.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_path, "prior-fit-diagnostics.csv"), row.names = FALSE)
write.csv(multiplier, file.path(output_path, "prior-multipliers.csv"), row.names = FALSE)
write.csv(ate, file.path(output_path, "prior-ates.csv"), row.names = FALSE)
write.csv(scales, file.path(output_path, "hierarchical-scales.csv"), row.names = FALSE)
write.csv(
  contrasts, file.path(output_path, "prior-multiplier-contrasts.csv"),
  row.names = FALSE
)

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
  file.path(output_path, "main-core-report-distance-prior-multipliers.tex")
)

contrast_columns <- c("500m", "1500m", "2500m", "500--2500m", "Minimum")
primary <- contrasts[
  contrasts$contrast == "No Signal - Any Signal",
]
contrast_lines <- c(
  "\\begin{tabular}{lccccc}", "\\toprule",
  "Specification & 500m & 1500m & 2500m & Finite & Grid minimum \\\\ ",
  "\\midrule"
)
for (specification_index in seq_len(nrow(specifications))) {
  cells <- primary[primary$id == specifications$id[specification_index], ]
  cells <- cells[match(contrast_columns, cells$distance), ]
  contrast_lines <- c(
    contrast_lines,
    paste0(
      specifications$label[specification_index], " & ",
      paste(fmt(cells$median), collapse = " & "), " \\\\"
    ),
    paste0(" & ", paste(ci(cells$lower, cells$upper), collapse = " & "),
           " \\\\ ")
  )
}
writeLines(
  c(contrast_lines, "\\bottomrule", "\\end{tabular}"),
  file.path(output_path, "main-core-report-distance-contrasts.tex")
)
print(diagnostics)
