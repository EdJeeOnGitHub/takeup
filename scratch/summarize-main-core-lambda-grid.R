#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
manifest_path <- main_core_option_value(args, "--manifest")
output_path <- main_core_option_value(
  args, "--output-path", "temp-data/main-core-lambda-identification"
)
if (is.null(manifest_path) || !file.exists(manifest_path)) {
  stop("--manifest must name an existing CSV.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(ggplot2)
})

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
treatment_labels <- c("Control", "Ink", "Calendar", "Bracelet")
distance_labels <- c("500m", "1500m", "2500m")
summarize_draw <- function(value) {
  c(
    mean = mean(value), sd = sd(value),
    q025 = unname(quantile(value, 0.025)),
    q50 = unname(quantile(value, 0.5)),
    q975 = unname(quantile(value, 0.975))
  )
}

all_estimands <- list()
all_diagnostics <- list()
for (row_index in seq_len(nrow(manifest))) {
  spec <- manifest[row_index, ]
  fit_dir <- file.path(output_path, "fits", spec$label)
  gq_dir <- file.path(output_path, "gq", spec$label)
  fit_csvs <- sort(list.files(
    fit_dir, pattern = "\\.csv$", full.names = TRUE
  ))
  gq_csvs <- sort(list.files(
    gq_dir, pattern = "\\.csv$", full.names = TRUE
  ))
  if (length(fit_csvs) != 4L || length(gq_csvs) != 4L) {
    warning("Skipping incomplete specification: ", spec$label)
    next
  }
  fit <- read_cmdstan_csv(fit_csvs)
  gq <- read_cmdstan_csv(gq_csvs)
  fit_summary <- summarise_draws(fit$post_warmup_draws)
  sampler_diag <- as_draws_df(fit$post_warmup_sampler_diagnostics)
  all_diagnostics[[length(all_diagnostics) + 1L]] <- data.frame(
    spec_id = spec$spec_id,
    label = spec$label,
    structure = spec$structure,
    prior_sd = spec$lambda_log_ratio_sd_prior,
    max_rhat = max(fit_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(fit_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(fit_summary$ess_tail, na.rm = TRUE),
    divergences = sum(sampler_diag$divergent__, na.rm = TRUE),
    max_treedepth = sum(sampler_diag$treedepth__ >= 12, na.rm = TRUE)
  )

  draws <- as_draws_df(gq$generated_quantities)
  add_estimand <- function(name, value, treatment = NA_character_,
                           distance = NA_character_) {
    summary <- summarize_draw(value)
    all_estimands[[length(all_estimands) + 1L]] <<- data.frame(
      spec_id = spec$spec_id, label = spec$label,
      structure = spec$structure,
      prior_sd = spec$lambda_log_ratio_sd_prior,
      estimand = name, treatment = treatment, distance = distance,
      mean = summary["mean"], sd = summary["sd"],
      q025 = summary["q025"], q50 = summary["q50"],
      q975 = summary["q975"]
    )
  }
  for (treatment_index in 1:4) {
    add_estimand(
      "lambda",
      draws[[sprintf("core_compact_signal_lambda[%d]", treatment_index)]],
      treatment_labels[treatment_index]
    )
    for (distance_index in 1:3) {
      add_estimand(
        "social_multiplier",
        -draws[[sprintf(
          "core_compact_sm_rescaled[%d,%d]",
          distance_index,
          treatment_index
        )]],
        treatment_labels[treatment_index],
        distance_labels[distance_index]
      )
    }
  }
  add_estimand(
    "signal_log_lambda_contrast",
    draws$core_compact_signal_vs_no_signal_log_lambda
  )
  for (distance_index in 1:3) {
    control <- draws[[sprintf(
      "core_compact_takeup_level[%d,1]", distance_index
    )]]
    for (treatment_index in 2:4) {
      add_estimand(
        "takeup_ate",
        draws[[sprintf(
          "core_compact_takeup_level[%d,%d]",
          distance_index,
          treatment_index
        )]] - control,
        treatment_labels[treatment_index],
        c("Combined", "Close", "Far")[distance_index]
      )
    }
  }
}

estimands <- bind_rows(all_estimands)
diagnostics <- bind_rows(all_diagnostics)
if (nrow(estimands) == 0L) stop("No complete specifications found.", call. = FALSE)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(
  estimands, file.path(output_path, "lambda-estimand-summary.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics, file.path(output_path, "lambda-fit-diagnostics.csv"),
  row.names = FALSE
)
contrast_contraction <- estimands |>
  filter(estimand == "signal_log_lambda_contrast") |>
  mutate(
    contrast_prior_sd = case_when(
      structure == "grouped" ~ prior_sd,
      structure == "arm" ~ prior_sd / sqrt(2),
      TRUE ~ NA_real_
    ),
    posterior_to_prior_sd = sd / contrast_prior_sd
  )
write.csv(
  contrast_contraction,
  file.path(output_path, "lambda-prior-contraction.csv"),
  row.names = FALSE
)

plot_data <- estimands |>
  filter(
    structure != "common",
    estimand == "signal_log_lambda_contrast" |
      (estimand == "social_multiplier" & distance == "1500m")
  ) |>
  mutate(
    series = ifelse(
      estimand == "signal_log_lambda_contrast",
      "log lambda: Any Signal - No Signal",
      paste0("Multiplier: ", treatment)
    ),
    structure = recode(
      structure,
      grouped = "Grouped lambda",
      arm = "Arm-specific lambda"
    )
  )
prior_plot <- ggplot(
  plot_data,
  aes(x = prior_sd, y = q50, ymin = q025, ymax = q975, color = structure)
) +
  geom_hline(yintercept = 0, color = "grey75", size = 0.3) +
  geom_errorbar(width = 0.025, position = position_dodge(width = 0.035)) +
  geom_point(position = position_dodge(width = 0.035)) +
  facet_wrap(~series, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(plot_data$prior_sd))) +
  labs(x = "Prior SD of pairwise log-lambda ratio", y = NULL, color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_path, "main-core-lambda-prior-sensitivity.pdf"),
  prior_plot,
  width = 7.2,
  height = 5.2
)
contraction_plot <- contrast_contraction |>
  filter(structure != "common") |>
  mutate(structure = recode(
    structure, grouped = "Grouped lambda", arm = "Arm-specific lambda"
  )) |>
  ggplot(aes(prior_sd, posterior_to_prior_sd, color = structure)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = sort(unique(plot_data$prior_sd))) +
  labs(
    x = "Prior SD of pairwise log-lambda ratio",
    y = "Posterior SD / prior SD", color = NULL
  ) + theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(
  file.path(output_path, "main-core-lambda-prior-contraction.pdf"),
  contraction_plot, width = 5.2, height = 3.5
)

fmt <- function(x) sprintf("%.2f", x)
ci <- function(lo, hi) paste0("(", fmt(lo), ", ", fmt(hi), ")")
selected <- estimands |>
  filter(
    structure == "common" |
      abs(prior_sd - 0.25) < 1e-10
  )
model_order <- c("common", "grouped", "arm")
lines <- c(
  "\\begin{tabular}{llccc}",
  "\\toprule",
  "Model & Treatment & 500m & 1500m & 2500m \\\\",
  "\\midrule"
)
for (structure_name in model_order) {
  block <- selected |>
    filter(
      structure == structure_name,
      estimand == "social_multiplier"
    )
  if (nrow(block) == 0L) next
  model_label <- c(
    common = "Common $\\lambda$",
    grouped = "Signal-group $\\lambda$",
    arm = "Arm-specific $\\lambda$"
  )[[structure_name]]
  for (treatment_name in treatment_labels) {
    cells <- block |>
      filter(treatment == treatment_name) |>
      arrange(match(distance, distance_labels))
    lines <- c(
      lines,
      paste0(
        if (treatment_name == treatment_labels[1]) model_label else "",
        " & ", treatment_name, " & ",
        paste(fmt(cells$q50), collapse = " & "), " \\\\"
      ),
      paste0(
        " &  & ",
        paste(ci(cells$q025, cells$q975), collapse = " & "), " \\\\"
      )
    )
  }
  if (structure_name != tail(model_order, 1)) lines <- c(lines, "\\midrule")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(
  lines,
  file.path(output_path, "main-core-lambda-nested-multipliers.tex")
)

ate <- selected |>
  filter(estimand == "takeup_ate")
lines <- c(
  "\\begin{tabular}{llccc}", "\\toprule",
  "Model & Treatment & Combined & Close & Far \\\\", "\\midrule"
)
for (structure_name in model_order) {
  block <- ate |> filter(structure == structure_name)
  if (nrow(block) == 0L) next
  model_label <- c(
    common = "Common $\\lambda$", grouped = "Signal-group $\\lambda$",
    arm = "Arm-specific $\\lambda$"
  )[[structure_name]]
  for (treatment_name in c("Bracelet", "Calendar", "Ink")) {
    cells <- block |>
      filter(treatment == treatment_name) |>
      arrange(match(distance, c("Combined", "Close", "Far")))
    lines <- c(
      lines,
      paste0(
        if (treatment_name == "Bracelet") model_label else "",
        " & ", treatment_name, " & ",
        paste(fmt(cells$q50), collapse = " & "), " \\\\"
      ),
      paste0(
        " &  & ",
        paste(ci(cells$q025, cells$q975), collapse = " & "), " \\\\"
      )
    )
  }
  if (structure_name != tail(model_order, 1)) lines <- c(lines, "\\midrule")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(output_path, "main-core-lambda-nested-ates.tex"))

prior_lines <- c(
  "\\begin{tabular}{llrrrrrr}", "\\toprule",
  "Model & Prior SD & $\\log(\\lambda_S/\\lambda_N)$ & Contraction & Control & Ink & Calendar & Bracelet \\\\",
  "\\midrule"
)
for (row_index in seq_len(nrow(contrast_contraction))) {
  contrast_row <- contrast_contraction[row_index, ]
  if (contrast_row$structure == "common") next
  sm <- estimands |>
    filter(
      label == contrast_row$label,
      estimand == "social_multiplier", distance == "1500m"
    ) |>
    arrange(match(treatment, treatment_labels))
  prior_lines <- c(prior_lines, paste0(
    ifelse(contrast_row$structure == "grouped", "Grouped", "Arm-specific"),
    " & ", sprintf("%.2f", contrast_row$prior_sd),
    " & ", fmt(contrast_row$q50),
    " & ", fmt(contrast_row$posterior_to_prior_sd),
    " & ", paste(fmt(sm$q50), collapse = " & "), " \\\\"
  ))
}
prior_lines <- c(prior_lines, "\\bottomrule", "\\end{tabular}")
writeLines(
  prior_lines,
  file.path(output_path, "main-core-lambda-prior-sensitivity.tex")
)

print(diagnostics)
