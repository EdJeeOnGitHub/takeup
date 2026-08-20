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
      "core_compact_sm_finite",
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

  finite_multiplier <- draws |>
    select(draw_id, starts_with("core_compact_sm_finite[")) |>
    pivot_longer(-draw_id, names_to = "variable", values_to = "value") |>
    extract(variable, "treatment", "\\[(\\d+)\\]", convert = TRUE) |>
    transmute(
      draw_id, estimand = "finite_multiplier", subgroup = "500--2500m",
      treatment = treatment_labels[treatment], value
    )

  shock_sd <- draws |>
    transmute(
      draw_id, estimand = "cluster_shock_sd", subgroup = "All",
      treatment = "All", value = core_compact_cluster_shock_sd
    )

  bind_rows(
    takeup, incentive_ate, distance_ate, multiplier, finite_multiplier,
    shock_sd
  ) |>
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

contrast_draws <- all_draws |>
  filter(estimand %in% c("social_multiplier", "finite_multiplier")) |>
  select(specification, replicate, draw_id, estimand, subgroup, treatment, value) |>
  pivot_wider(names_from = treatment, values_from = value) |>
  filter(if_all(all_of(treatment_labels), is.finite)) |>
  mutate(
    `No Signal - Any Signal` = (Control + Calendar - Ink - Bracelet) / 2,
    `Control - Bracelet` = Control - Bracelet,
    `Calendar - Bracelet` = Calendar - Bracelet,
    `Control - Ink` = Control - Ink
  ) |>
  pivot_longer(
    all_of(c(
      "No Signal - Any Signal", "Control - Bracelet",
      "Calendar - Bracelet", "Control - Ink"
    )),
    names_to = "contrast", values_to = "value"
  ) |>
  select(
    specification, replicate, draw_id, estimand, subgroup, contrast, value
  )

grid_minimum <- contrast_draws |>
  filter(
    estimand == "social_multiplier",
    contrast == "No Signal - Any Signal",
    subgroup %in% distance_labels
  ) |>
  group_by(specification, replicate, draw_id, contrast) |>
  filter(n_distinct(subgroup) == length(distance_labels)) |>
  summarise(
    estimand = "grid_minimum", subgroup = "Minimum",
    value = min(value), .groups = "drop"
  )
contrast_draws <- bind_rows(contrast_draws, grid_minimum)
contrast_summary <- contrast_draws |>
  group_by(specification, estimand, subgroup, contrast) |>
  summarise(
    estimate = median(value),
    conf_low = quantile(value, 0.025),
    conf_high = quantile(value, 0.975),
    share_positive = mean(value > 0),
    draws_or_replications = n_distinct(replicate),
    .groups = "drop"
  )
write.csv(
  contrast_draws, file.path(output_path, "multiplier-contrast-draws.csv"),
  row.names = FALSE
)
write.csv(
  contrast_summary,
  file.path(output_path, "multiplier-contrast-summary.csv"), row.names = FALSE
)

if (length(status_files) > 0L) {
  status_detail <- map_dfr(
    status_files, read.csv, stringsAsFactors = FALSE
  )
  valid_weight <- !is.na(status_detail$weight_file) &
    file.exists(status_detail$weight_file)
  status_detail$weight_hash <- NA_character_
  status_detail$weight_hash[valid_weight] <- unname(tools::md5sum(
    status_detail$weight_file[valid_weight]
  ))
  status_audit <- status_detail |>
    group_by(method) |>
    summarise(
      attempted = n_distinct(replicate),
      successful = n_distinct(replicate[status == "complete"]),
      failed = n_distinct(replicate[status != "complete"]),
      included = n_distinct(replicate[
        status == "complete" & file.exists(gq_csv)
      ]),
      unique_weight_vectors = n_distinct(weight_hash[!is.na(weight_hash)]),
      duplicate_weight_vectors =
        sum(duplicated(weight_hash[!is.na(weight_hash)])),
      .groups = "drop"
    )
  write.csv(
    status_audit, file.path(output_path, "bootstrap-run-audit.csv"),
    row.names = FALSE
  )
}

multiplier_summary <- summary |>
  filter(estimand == "social_multiplier") |>
  mutate(
    specification = recode(
      specification,
      "Exponential cluster weights" =
        "Cluster weighted-likelihood bootstrap"
    ),
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
  "Baseline structural model", "Cluster weighted-likelihood bootstrap",
  "Cluster random shock", "Cluster random shock: prior SD 0.25"
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

primary_contrast <- contrast_summary |>
  filter(contrast == "No Signal - Any Signal") |>
  mutate(
    specification = recode(
      specification,
      "Exponential cluster weights" =
        "Cluster weighted-likelihood bootstrap"
    ),
    subgroup = factor(
      subgroup,
      levels = c(distance_labels, "500--2500m", "Minimum")
    ),
    specification = factor(specification, levels = spec_order)
  ) |>
  arrange(specification, subgroup)
contrast_path <- file.path(output_path, "main-core-multiplier-contrasts.tex")
contrast_con <- file(contrast_path, open = "wt")
on.exit(close(contrast_con), add = TRUE)
writeLines(c(
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Specification & 500m & 1500m & 2500m & Finite & Grid minimum \\\\ ",
  "\\midrule"
), contrast_con)
for (spec in available_specs) {
  row <- primary_contrast |>
    filter(specification == spec) |>
    arrange(subgroup)
  if (nrow(row) != 5L) next
  writeLines(sprintf(
    "%s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\ ",
    escape_tex(spec), row$estimate[1], row$estimate[2], row$estimate[3],
    row$estimate[4], row$estimate[5]
  ), contrast_con)
  writeLines(sprintf(
    " & [%.2f, %.2f] & [%.2f, %.2f] & [%.2f, %.2f] & [%.2f, %.2f] & [%.2f, %.2f] \\\\ ",
    row$conf_low[1], row$conf_high[1], row$conf_low[2], row$conf_high[2],
    row$conf_low[3], row$conf_high[3], row$conf_low[4], row$conf_high[4],
    row$conf_low[5], row$conf_high[5]
  ), contrast_con)
}
writeLines(c("\\bottomrule", "\\end{tabular}"), contrast_con)
message("Wrote summaries, table, and figure under ", output_path)
