#!/usr/bin/env Rscript

library(tidyverse)
library(fixest)

script_options <- docopt::docopt(
  "Usage:
  prior-deworming-robustness-table.R [options]

Options:
  --table-output-path=<path>  Table output path [default: presentations/rf-tables/main-specs]
",
  args = if (interactive()) "" else commandArgs(trailingOnly = TRUE)
)

source(file.path("scratch", "reduced-form-setup.R"))

table_output_path <- script_options$table_output_path
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)

sd_of_dist <- sd(analysis_data$cluster.dist.to.pot, na.rm = TRUE)

cluster_expected_dist <- readr::read_csv(
  file.path("data", "cluster_expected_dist.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    cluster.id = as.character(cluster.id),
    mu_d = dist / sd_of_dist
  )

baseline_prior_deworming <- baseline_worm %>%
  mutate(cluster.id = as.character(cluster.id)) %>%
  group_by(cluster.id) %>%
  summarise(
    baseline_ever_dewormed = mean(treated_lgl, na.rm = TRUE),
    baseline_past_year_dewormed = mean(treated_past_year, na.rm = TRUE),
    baseline_n = n(),
    .groups = "drop"
  )

analysis_for_prior <- analysis_data %>%
  mutate(
    cluster.id = as.character(cluster.id),
    assigned_dist_group = factor(assigned_dist_group, levels = c("close", "far")),
    public_signal = if_else(
      assigned_treatment %in% c("bracelet", "ink"),
      "Public signal",
      "No public signal"
    ),
    public_signal = factor(public_signal, levels = c("No public signal", "Public signal")),
    dewormed_monitored = as.numeric(dewormed),
    dewormed_self_reported = as.numeric(dewormed.reported)
  ) %>%
  left_join(baseline_prior_deworming, by = "cluster.id") %>%
  left_join(cluster_expected_dist, by = "cluster.id") %>%
  filter(
    !is.na(baseline_ever_dewormed),
    !is.na(baseline_past_year_dewormed),
    !is.na(mu_d)
  )

fit_signal_model <- function(data, outcome, prior_controls = FALSE) {
  rhs_terms <- c(
    "0 + public_signal * assigned_dist_group",
    "female",
    "age.census",
    "mu_d",
    if (prior_controls) c("baseline_ever_dewormed", "baseline_past_year_dewormed")
  )
  rhs <- paste(rhs_terms, collapse = " + ")
  feols(
    as.formula(paste(outcome, "~", rhs, "| county")),
    data = data,
    cluster = ~ cluster.id
  )
}

fits <- list(
  "Monitored" = fit_signal_model(analysis_for_prior, "dewormed_monitored", FALSE),
  "Monitored + prior" = fit_signal_model(analysis_for_prior, "dewormed_monitored", TRUE),
  "Self-reported" = fit_signal_model(analysis_for_prior, "dewormed_self_reported", FALSE),
  "Self-reported + prior" = fit_signal_model(analysis_for_prior, "dewormed_self_reported", TRUE)
)

fit_metadata <- tibble(
  model = names(fits),
  outcome = c(
    "dewormed_monitored",
    "dewormed_monitored",
    "dewormed_self_reported",
    "dewormed_self_reported"
  )
) %>%
  mutate(
    observations = map_int(outcome, ~ sum(!is.na(analysis_for_prior[[.x]]))),
    clusters = map_int(
      outcome,
      ~ n_distinct(analysis_for_prior$cluster.id[!is.na(analysis_for_prior[[.x]])])
    )
  )

term_labels <- c(
  "public_signalPublic signal:assigned_dist_groupfar" = "Public signal $\\times$ Far",
  "assigned_dist_groupfar" = "Far",
  "baseline_ever_dewormed" = "Baseline ever dewormed share",
  "baseline_past_year_dewormed" = "Baseline past-year dewormed share"
)

format_estimate <- function(fit, term) {
  coefs <- coef(fit)
  if (!term %in% names(coefs)) return("")
  se <- se(fit)
  sprintf("%.3f\n(%.3f)", unname(coefs[term]), unname(se[term]))
}

format_n <- function(fit) format(nobs(fit), big.mark = ",", scientific = FALSE)

format_count <- function(x) {
  format(
    x,
    big.mark = ",",
    scientific = FALSE
  )
}

make_body_row <- function(label, values) {
  paste0(label, " & ", paste(values, collapse = " & "), " \\\\")
}

estimate_rows <- imap(term_labels, function(label, term) {
  make_body_row(
    label,
    map_chr(fits, ~ format_estimate(.x, term) %>% str_replace("\n", " "))
  )
}) %>%
  unlist(use.names = FALSE)

metadata_rows <- c(
  "\\midrule",
  make_body_row("Observations", map_chr(fit_metadata$observations, format_count)),
  make_body_row("Clusters", map_chr(fit_metadata$clusters, format_count)),
  make_body_row("County fixed effects", rep("Yes", length(fits))),
  make_body_row("Main covariates", rep("Yes", length(fits))),
  make_body_row("Baseline prior-deworming shares", c("No", "Yes", "No", "Yes"))
)

table_lines <- c(
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  paste0(" & ", paste0("(", seq_along(fits), ") ", names(fits), collapse = " & "), " \\\\"),
  "\\midrule",
  unname(estimate_rows),
  metadata_rows,
  "\\bottomrule",
  "\\end{tabular}"
)

writeLines(
  table_lines,
  file.path(table_output_path, "prior-deworming-robustness.tex")
)

message("Wrote ", file.path(table_output_path, "prior-deworming-robustness.tex"))
