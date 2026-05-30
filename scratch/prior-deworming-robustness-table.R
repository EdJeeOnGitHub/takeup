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
    assigned_treatment = factor(
      assigned_treatment,
      levels = c("control", "calendar", "bracelet", "ink")
    ),
    assigned_dist_group = factor(assigned_dist_group, levels = c("close", "far")),
    dewormed_monitored = as.numeric(dewormed)
  ) %>%
  left_join(baseline_prior_deworming, by = "cluster.id") %>%
  left_join(cluster_expected_dist, by = "cluster.id") %>%
  filter(
    !is.na(baseline_ever_dewormed),
    !is.na(baseline_past_year_dewormed),
    !is.na(mu_d)
  )

prior_terms <- c("baseline_ever_dewormed", "baseline_past_year_dewormed")
base_controls <- c("female", "age.census", "mu_d", prior_terms)

fit_arm_model <- function(data, outcome, distance_group = NULL) {
  model_data <- data %>% filter(!is.na(.data[[outcome]]))
  if (!is.null(distance_group)) {
    model_data <- model_data %>% filter(assigned_dist_group == distance_group)
  }

  rhs <- paste(c("assigned_treatment", base_controls), collapse = " + ")
  feols(
    as.formula(paste(outcome, "~", rhs, "| county")),
    data = model_data,
    cluster = ~ cluster.id
  )
}

fit_distance_model <- function(data, outcome) {
  model_data <- data %>% filter(!is.na(.data[[outcome]]))
  rhs <- paste(c("assigned_treatment * assigned_dist_group", base_controls), collapse = " + ")
  feols(
    as.formula(paste(outcome, "~", rhs, "| county")),
    data = model_data,
    cluster = ~ cluster.id
  )
}

fit_bundle <- function(data, outcome) {
  list(
    combined = fit_arm_model(data, outcome),
    close = fit_arm_model(data, outcome, "close"),
    far = fit_arm_model(data, outcome, "far"),
    distance = fit_distance_model(data, outcome)
  )
}

fits <- list(
  "Monitored take-up" = fit_bundle(analysis_for_prior, "dewormed_monitored")
)

format_count <- function(x) {
  format(
    x,
    big.mark = ",",
    scientific = FALSE
  )
}

stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
}

format_estimate <- function(fit, term) {
  ct <- coeftable(fit)
  if (!term %in% rownames(ct)) return("")
  estimate <- ct[term, "Estimate"]
  std_error <- ct[term, "Std. Error"]
  p_value <- ct[term, "Pr(>|t|)"]
  sprintf(
    "\\makecell[c]{%.3f%s\\\\{(%.3f)}}",
    estimate,
    stars(p_value),
    std_error
  )
}

control_summary <- function(data, outcome, distance_group = NULL) {
  model_data <- data %>%
    filter(assigned_treatment == "control", !is.na(.data[[outcome]]))
  if (!is.null(distance_group)) {
    model_data <- model_data %>% filter(assigned_dist_group == distance_group)
  }

  fit <- feols(
    as.formula(paste(outcome, "~ 1")),
    data = model_data,
    cluster = ~ cluster.id
  )

  sprintf(
    "\\makecell[c]{%.3f\\\\{(%.3f)}}",
    unname(coef(fit)[["(Intercept)"]]),
    unname(se(fit)[["(Intercept)"]])
  )
}

term_for_arm <- function(arm, column) {
  if (arm == "control") {
    if (column == "distance") return("assigned_dist_groupfar")
    return(NA_character_)
  }

  treatment_term <- paste0("assigned_treatment", arm)
  if (column == "distance") {
    paste0(treatment_term, ":assigned_dist_groupfar")
  } else {
    treatment_term
  }
}

cell_for_arm <- function(fit_set, data, outcome, arm, column) {
  if (arm == "control" && column != "distance") {
    return(control_summary(
      data,
      outcome,
      distance_group = if (column == "combined") NULL else column
    ))
  }

  term <- term_for_arm(arm, column)
  format_estimate(fit_set[[column]], term)
}

make_body_row <- function(label, values) {
  paste0(label, " & ", paste(values, collapse = " & "), " \\\\")
}

make_panel <- function(panel_name, fit_set, outcome) {
  columns <- c("combined", "close", "far", "distance")
  arm_labels <- c(
    bracelet = "Bracelet",
    calendar = "Calendar",
    ink = "Ink",
    control = "Control"
  )

  treatment_rows <- imap(arm_labels, function(label, arm) {
    make_body_row(
      label,
      map_chr(columns, ~ cell_for_arm(fit_set, analysis_for_prior, outcome, arm, .x))
    )
  }) %>%
    unlist(use.names = FALSE)

  covariate_rows <- c(
    "\\addlinespace[0.3em]",
    make_body_row(
      "Baseline ever dewormed share",
      c(
        format_estimate(fit_set$combined, "baseline_ever_dewormed"),
        format_estimate(fit_set$close, "baseline_ever_dewormed"),
        format_estimate(fit_set$far, "baseline_ever_dewormed"),
        format_estimate(fit_set$distance, "baseline_ever_dewormed")
      )
    ),
    make_body_row(
      "Baseline past-year dewormed share",
      c(
        format_estimate(fit_set$combined, "baseline_past_year_dewormed"),
        format_estimate(fit_set$close, "baseline_past_year_dewormed"),
        format_estimate(fit_set$far, "baseline_past_year_dewormed"),
        format_estimate(fit_set$distance, "baseline_past_year_dewormed")
      )
    )
  )

  metadata_rows <- c(
    "\\addlinespace[0.3em]",
    make_body_row("Observations", map_chr(fit_set, ~ format_count(nobs(.x)))),
    make_body_row(
      "Clusters",
      map_chr(
        list(NULL, "close", "far", NULL),
        function(distance_group) {
          model_data <- analysis_for_prior %>% filter(!is.na(.data[[outcome]]))
          if (!is.null(distance_group)) {
            model_data <- model_data %>% filter(assigned_dist_group == distance_group)
          }
          format_count(n_distinct(model_data$cluster.id))
        }
      )
    ),
    make_body_row("County fixed effects", rep("Yes", length(columns))),
    make_body_row("Main covariates", rep("Yes", length(columns))),
    make_body_row("Baseline prior-deworming covariates", rep("Yes", length(columns)))
  )

  c(
    paste0("\\multicolumn{5}{l}{\\textbf{", panel_name, "}} \\\\"),
    "\\addlinespace[0.2em]",
    treatment_rows,
    covariate_rows,
    metadata_rows
  )
}

table_lines <- c(
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & Combined & Close & Far & Far--Close \\\\",
  make_panel("Monitored take-up", fits[["Monitored take-up"]], "dewormed_monitored"),
  "\\bottomrule",
  "\\end{tabular}"
)

writeLines(
  table_lines,
  file.path(table_output_path, "prior-deworming-robustness.tex")
)

message("Wrote ", file.path(table_output_path, "prior-deworming-robustness.tex"))
