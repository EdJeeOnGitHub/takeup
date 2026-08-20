#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(ggplot2)
})

baseline_gq_dir <- option_value("--baseline-gq-dir")
mixture_fit_dir <- option_value("--mixture-fit-dir")
mixture_gq_dir <- option_value("--mixture-gq-dir")
output_path <- option_value(
  "--output-path", "temp-data/main-core-finite-mixture-robustness"
)
required_dirs <- c(baseline_gq_dir, mixture_fit_dir, mixture_gq_dir)
if (any(vapply(
  required_dirs, function(path) is.null(path) || !dir.exists(path), logical(1)
))) {
  stop("Existing baseline GQ, mixture fit, and mixture GQ directories are required.",
       call. = FALSE)
}

csvs <- function(path, pattern = "\\.csv$") {
  sort(list.files(path, pattern = pattern, full.names = TRUE))
}
baseline_gq_csvs <- csvs(baseline_gq_dir)
mixture_fit_csvs <- csvs(mixture_fit_dir, "finite-mixture-chain[1-4]-.*\\.csv$")
mixture_gq_csvs <- csvs(mixture_gq_dir)
if (length(baseline_gq_csvs) != 4L || length(mixture_fit_csvs) != 4L ||
    length(mixture_gq_csvs) != 4L) {
  stop("Expected four CSVs for each fit/GQ input.", call. = FALSE)
}

baseline <- as_draws_df(read_cmdstan_csv(baseline_gq_csvs)$generated_quantities)
mixture_gq <- as_draws_df(read_cmdstan_csv(mixture_gq_csvs)$generated_quantities)
mixture_fit <- read_cmdstan_csv(mixture_fit_csvs)
fit_draws <- as_draws_df(mixture_fit$post_warmup_draws)
fit_summary <- summarise_draws(mixture_fit$post_warmup_draws)
sampler <- as_draws_df(mixture_fit$post_warmup_sampler_diagnostics)

treatments <- c("Control", "Ink", "Calendar", "Bracelet")
distances <- c("500m", "1500m", "2500m")
summarize_value <- function(value) c(
  median = median(value),
  q025 = unname(quantile(value, 0.025)),
  q975 = unname(quantile(value, 0.975))
)

rows <- list()
for (specification in c("Gaussian", "Two-component mixture")) {
  draws <- if (specification == "Gaussian") baseline else mixture_gq
  for (treatment_index in seq_along(treatments)) {
    for (distance_index in seq_along(distances)) {
      value <- -draws[[sprintf(
        "core_compact_sm_rescaled[%d,%d]", distance_index, treatment_index
      )]]
      summary <- summarize_value(value)
      rows[[length(rows) + 1L]] <- data.frame(
        specification = specification,
        treatment = treatments[treatment_index],
        distance = distances[distance_index],
        median = summary[["median"]], q025 = summary[["q025"]],
        q975 = summary[["q975"]]
      )
    }
  }
}
results <- bind_rows(rows)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output_path, "finite-mixture-multiplier-summary.csv"),
          row.names = FALSE)

mixture_weight <- fit_draws[["core_finite_mixture_weight[1]"]]
between_share <- fit_draws[["core_finite_mixture_between_share[1]"]]
separation <- sqrt(between_share / (mixture_weight * (1 - mixture_weight)))
shape_draws <- data.frame(
  mixture_weight = mixture_weight,
  between_variance_share = between_share,
  lower_component_mean = -(1 - mixture_weight) * separation,
  upper_component_mean = mixture_weight * separation,
  within_component_sd = sqrt(1 - between_share)
)
shape_summary <- bind_rows(lapply(names(shape_draws), function(variable) {
  summary <- summarize_value(shape_draws[[variable]])
  data.frame(
    variable = variable, median = summary[["median"]],
    q025 = summary[["q025"]], q975 = summary[["q975"]]
  )
}))
write.csv(shape_summary, file.path(output_path, "finite-mixture-shape-summary.csv"),
          row.names = FALSE)

diagnostics <- data.frame(
  max_rhat_variable = fit_summary$variable[which.max(fit_summary$rhat)],
  max_rhat = max(fit_summary$rhat, na.rm = TRUE),
  min_ess_bulk_variable = fit_summary$variable[which.min(fit_summary$ess_bulk)],
  min_ess_bulk = min(fit_summary$ess_bulk, na.rm = TRUE),
  min_ess_tail_variable = fit_summary$variable[which.min(fit_summary$ess_tail)],
  min_ess_tail = min(fit_summary$ess_tail, na.rm = TRUE),
  divergences = sum(sampler$divergent__, na.rm = TRUE),
  max_treedepth = sum(sampler$treedepth__ >= 12, na.rm = TRUE)
)
write.csv(diagnostics, file.path(output_path, "finite-mixture-fit-diagnostics.csv"),
          row.names = FALSE)

# Plot the estimated standardized intrinsic-motivation distribution. The
# parameterization fixes its mean at zero and variance at one; the mixture
# flexibly changes skewness, shoulders, and modality relative to N(0,1).
density_grid <- seq(-3.5, 3.5, length.out = 281L)
density_draws <- vapply(seq_len(nrow(shape_draws)), function(draw_index) {
  draw <- shape_draws[draw_index, ]
  draw$mixture_weight * dnorm(
    density_grid, draw$lower_component_mean, draw$within_component_sd
  ) + (1 - draw$mixture_weight) * dnorm(
    density_grid, draw$upper_component_mean, draw$within_component_sd
  )
}, numeric(length(density_grid)))
density_summary <- data.frame(
  intrinsic_motivation = density_grid,
  mixture_median = apply(density_draws, 1, median),
  mixture_q025 = apply(density_draws, 1, quantile, probs = 0.025),
  mixture_q975 = apply(density_draws, 1, quantile, probs = 0.975),
  gaussian_density = dnorm(density_grid)
)
median_shape <- vapply(shape_draws, median, numeric(1))
density_summary$lower_component <-
  median_shape[["mixture_weight"]] * dnorm(
    density_grid,
    median_shape[["lower_component_mean"]],
    median_shape[["within_component_sd"]]
  )
density_summary$upper_component <-
  (1 - median_shape[["mixture_weight"]]) * dnorm(
    density_grid,
    median_shape[["upper_component_mean"]],
    median_shape[["within_component_sd"]]
  )
write.csv(
  density_summary,
  file.path(output_path, "finite-mixture-density-summary.csv"),
  row.names = FALSE
)

density_plot <- ggplot(
  density_summary, aes(x = .data$intrinsic_motivation)
) +
  geom_ribbon(
    aes(ymin = .data$mixture_q025, ymax = .data$mixture_q975),
    fill = "#8172B2", alpha = 0.20
  ) +
  geom_line(
    aes(y = .data$mixture_median, color = "Estimated mixture"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = .data$gaussian_density, color = "Standard normal"),
    linewidth = 0.9, linetype = "longdash"
  ) +
  geom_line(
    aes(y = .data$lower_component, color = "Weighted components"),
    linewidth = 0.55, linetype = "dotted"
  ) +
  geom_line(
    aes(y = .data$upper_component, color = "Weighted components"),
    linewidth = 0.55, linetype = "dotted"
  ) +
  scale_color_manual(values = c(
    "Estimated mixture" = "#6A51A3",
    "Standard normal" = "#333333",
    "Weighted components" = "#4C72B0"
  )) +
  labs(
    x = "Standardized intrinsic motivation",
    y = "Density",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
ggsave(
  file.path(output_path, "finite-mixture-density.pdf"),
  density_plot, width = 6.4, height = 4.2
)

fmt <- function(x) sprintf("%.2f", x)
ci <- function(lo, hi) paste0("(", fmt(lo), ", ", fmt(hi), ")")
lines <- c(
  "\\begin{tabular}{lccc}", "\\toprule",
  "Treatment & 500m & 1500m & 2500m \\\\",
  "\\midrule"
)
for (specification_name in c("Gaussian", "Two-component mixture")) {
  block <- results[
    results[["specification"]] == specification_name, , drop = FALSE
  ]
  panel <- if (specification_name == "Gaussian") {
    "Panel A: Gaussian intrinsic motivation"
  } else {
    "Panel B: Two-component Gaussian-mixture intrinsic motivation"
  }
  lines <- c(lines, paste0("\\multicolumn{4}{l}{\\textit{", panel, "}} \\\\"))
  for (treatment_name in c("Bracelet", "Calendar", "Ink", "Control")) {
    cells <- block |>
      filter(.data$treatment == .env$treatment_name) |>
      arrange(match(.data$distance, distances))
    lines <- c(
      lines,
      paste0(treatment_name, " & ", paste(fmt(cells$median), collapse = " & "),
             " \\\\"),
      paste0(" & ", paste(ci(cells$q025, cells$q975), collapse = " & "),
             " \\\\ ")
    )
  }
  if (specification_name == "Gaussian") lines <- c(lines, "\\addlinespace")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(output_path, "main-core-finite-mixture-multipliers.tex"))
print(diagnostics)
print(shape_summary)
