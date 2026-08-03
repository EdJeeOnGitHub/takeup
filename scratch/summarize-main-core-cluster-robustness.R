#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(posterior)
})

baseline_csv <- option_value("--baseline-gq")
shock_csv <- option_value("--shock-gq")
shock_sensitivity_csv <- option_value("--shock-sensitivity-gq")
weighted_path <- option_value(
  "--weighted-path", "temp/main-core-weighted-modes"
)
output_path <- option_value(
  "--output-path", "temp-data/main-core-cluster-robustness"
)
table_path <- option_value(
  "--table-path",
  "presentations/tables/fit105/main-core-cluster-robustness.tex"
)
figure_path <- option_value(
  "--figure-path",
  "presentations/misc-figures/main-core-cluster-robustness.pdf"
)

split_csvs <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  paths <- strsplit(x, ",", fixed = TRUE)[[1L]]
  if (any(!file.exists(paths))) stop("Missing GQ CSV.", call. = FALSE)
  paths
}

treatment_labels <- c("Control", "Ink", "Calendar", "Bracelet")
distance_labels <- c("500m", "1500m", "2500m")
takeup_labels <- c("Combined", "Close", "Far")

read_compact <- function(files, specification, replicate = NULL) {
  if (length(files) == 0L) return(tibble())
  draws <- read_cmdstan_csv(
    files,
    variables = c(
      "core_compact_takeup_level",
      "core_compact_sm_rescaled",
      "core_compact_cluster_shock_sd"
    )
  )$generated_quantities |>
    as_draws_df() |>
    as.data.frame() |>
    mutate(draw_id = row_number())

  takeup <- draws |>
    select(draw_id, starts_with("core_compact_takeup_level[")) |>
    pivot_longer(-draw_id, names_to = "variable", values_to = "value") |>
    extract(variable, c("row", "treatment"), "\\[(\\d+),(\\d+)\\]", convert = TRUE) |>
    mutate(
      estimand = "takeup_level",
      subgroup = takeup_labels[row],
      treatment = treatment_labels[treatment]
    ) |>
    select(draw_id, estimand, subgroup, treatment, value)

  incentive_ate <- takeup |>
    select(draw_id, subgroup, treatment, value) |>
    left_join(
      takeup |> filter(treatment == "Control") |>
        select(draw_id, subgroup, control = value),
      by = c("draw_id", "subgroup")
    ) |>
    filter(treatment != "Control") |>
    transmute(
      draw_id, estimand = "incentive_ate", subgroup, treatment,
      value = value - control
    )

  distance_ate <- takeup |>
    filter(subgroup %in% c("Close", "Far")) |>
    select(draw_id, subgroup, treatment, value) |>
    pivot_wider(names_from = subgroup, values_from = value) |>
    transmute(
      draw_id, estimand = "far_minus_close", subgroup = "Far - Close",
      treatment, value = Far - Close
    )

  multiplier <- draws |>
    select(draw_id, starts_with("core_compact_sm_rescaled[")) |>
    pivot_longer(-draw_id, names_to = "variable", values_to = "value") |>
    extract(variable, c("row", "treatment"), "\\[(\\d+),(\\d+)\\]", convert = TRUE) |>
    transmute(
      draw_id, estimand = "social_multiplier",
      subgroup = distance_labels[row],
      treatment = treatment_labels[treatment],
      # Existing paper convention is the sign-reversed rescaled multiplier.
      value = -value
    )

  shock_sd <- draws |>
    transmute(
      draw_id, estimand = "cluster_shock_sd", subgroup = "All",
      treatment = "All", value = core_compact_cluster_shock_sd
    )

  bind_rows(takeup, incentive_ate, distance_ate, multiplier, shock_sd) |>
    mutate(
      specification = specification,
      replicate = if (is.null(replicate)) draw_id else replicate
    )
}

all_draws <- bind_rows(
  read_compact(split_csvs(baseline_csv), "Baseline structural model"),
  read_compact(split_csvs(shock_csv), "Cluster random shock"),
  read_compact(
    split_csvs(shock_sensitivity_csv),
    "Cluster random shock: prior SD 0.25"
  )
)

status_files <- Sys.glob(file.path(weighted_path, "*", "status.csv"))
if (length(status_files) > 0L) {
  completed_status <- map_dfr(
    status_files, read.csv, stringsAsFactors = FALSE
  ) |>
    filter(
      status == "complete",
      method %in% c("multinomial", "exponential"),
      file.exists(gq_csv)
    ) |>
    group_by(method) |>
    arrange(replicate, .by_group = TRUE) |>
    slice_head(n = 200L) |>
    ungroup()
  weighted_draws <- map_dfr(seq_len(nrow(completed_status)), function(i) {
    status <- completed_status[i, ]
    label <- if (status$method == "multinomial") {
      "Cluster bootstrap"
    } else if (status$method == "exponential") {
      "Exponential cluster weights"
    } else {
      "Baseline structural model"
    }
    read_compact(status$gq_csv, label, status$replicate)
  })
  all_draws <- bind_rows(all_draws, weighted_draws)
}
if (nrow(all_draws) == 0L) stop("No compact GQ inputs found.", call. = FALSE)

summary <- all_draws |>
  group_by(specification, estimand, subgroup, treatment) |>
  summarise(
    estimate = median(value, na.rm = TRUE),
    conf_low = quantile(value, 0.025, na.rm = TRUE),
    conf_high = quantile(value, 0.975, na.rm = TRUE),
    draws_or_replications = n_distinct(replicate),
    finite_draws_or_replications = n_distinct(replicate[is.finite(value)]),
    .groups = "drop"
  )

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
write.csv(all_draws, file.path(output_path, "estimand-draws.csv"), row.names = FALSE)
write.csv(summary, file.path(output_path, "estimand-summary.csv"), row.names = FALSE)

multiplier_summary <- summary |>
  filter(estimand == "social_multiplier") |>
  mutate(
    subgroup = factor(subgroup, levels = distance_labels),
    treatment = factor(treatment, levels = c("Control", "Ink", "Bracelet", "Calendar"))
  )

plot <- ggplot(
  multiplier_summary,
  aes(x = subgroup, y = estimate, ymin = conf_low, ymax = conf_high,
      color = specification, group = specification)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_wrap(~ treatment, ncol = 2) +
  labs(x = "Distance", y = "Social multiplier", color = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(figure_path, plot, width = 7.2, height = 5.2)

escape_tex <- function(x) gsub("_", "\\\\_", x, fixed = TRUE)
spec_order <- c(
  "Baseline structural model", "Cluster bootstrap", "Cluster random shock",
  "Cluster random shock: prior SD 0.25",
  "Exponential cluster weights"
)
table_data <- multiplier_summary |>
  mutate(specification = factor(specification, levels = spec_order)) |>
  arrange(specification, treatment, subgroup)

con <- file(table_path, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(c(
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Treatment & 500m & 1500m & 2500m \\\\",
  "\\midrule"
), con)
available_specs <- spec_order[spec_order %in% as.character(table_data$specification)]
for (spec in available_specs) {
  writeLines(sprintf(
    "\\multicolumn{4}{l}{\\textit{%s}} \\\\", escape_tex(as.character(spec))
  ), con)
  panel <- table_data |> filter(specification == spec)
  for (treat in levels(panel$treatment)) {
    row <- panel |> filter(treatment == treat) |> arrange(subgroup)
    if (nrow(row) != 3L) next
    writeLines(sprintf(
      "%s & %.2f & %.2f & %.2f \\\\",
      treat, row$estimate[1], row$estimate[2], row$estimate[3]
    ), con)
    writeLines(sprintf(
      " & [%.2f, %.2f] & [%.2f, %.2f] & [%.2f, %.2f] \\\\",
      row$conf_low[1], row$conf_high[1], row$conf_low[2], row$conf_high[2],
      row$conf_low[3], row$conf_high[3]
    ), con)
  }
  writeLines("\\addlinespace", con)
}
writeLines(c("\\bottomrule", "\\end{tabular}"), con)
message("Wrote summaries, table, and figure under ", output_path)
