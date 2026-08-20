#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(prefix, default) {
  hit <- args[str_starts(args, paste0(prefix, "="))]
  if (length(hit) == 0) return(default)
  str_sub(hit[[1]], str_length(prefix) + 2)
}

B <- as.integer(get_arg("--bootstrap-draws", "1000"))
bootstrap_seed <- as.integer(get_arg("--seed", "20260807"))
output_dir <- get_arg("--output-dir", "ref-reports/campaign-day")
appendix_dir <- get_arg("--appendix-dir", "appendix/structural-robustness")
if (is.na(B) || B < 1) stop("--bootstrap-draws must be a positive integer.")
if (is.na(bootstrap_seed)) stop("--seed must be an integer.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

clean_cluster_id <- function(x, label = "cluster ID") {
  x_chr <- str_trim(as.character(x))
  bad <- is.na(x_chr) | x_chr == "" | !str_detect(x_chr, "^[0-9]+$")
  if (any(bad)) stop(label, " contains missing or non-integer values.", call. = FALSE)
  as.integer(x_chr)
}

message("Reading the original-assignment analysis sample")
analysis_data <- read_csv(
  "temp-data/analysis-cluster-recentered-covariate-data.csv",
  show_col_types = FALSE,
  progress = FALSE
) %>%
  transmute(
    KEY.individ,
    original_cluster_id = clean_cluster_id(cluster.id.x, "analysis cluster ID"),
    cluster_id,
    county = factor(county),
    assigned_treatment = factor(
      assigned.treatment,
      levels = c("control", "ink", "calendar", "bracelet")
    ),
    dewormed = as.logical(dewormed),
    dewormed_day = as.integer(dewormed.day),
    female,
    age.census,
    mu_d
  )

original_assignment <- read_rds("data/rct_targetable_schools_2.0.rds") %>%
  as_tibble() %>%
  transmute(
    original_cluster_id = clean_cluster_id(pot.cluster.id, "assignment cluster ID"),
    assigned_dist_group = factor(assigned.dist.cat, levels = c("close", "far"))
  ) %>%
  distinct()

assert_true(
  nrow(original_assignment) == n_distinct(original_assignment$original_cluster_id),
  "Original distance assignment is not unique by cluster."
)

analysis_data <- analysis_data %>%
  left_join(original_assignment, by = "original_cluster_id", relationship = "many-to-one")

assert_true(nrow(analysis_data) == 9805, "Expected 9,805 observations.")
assert_true(!anyDuplicated(analysis_data$KEY.individ), "Individual identifiers are not unique.")
assert_true(n_distinct(analysis_data$original_cluster_id) == 144, "Expected 144 clusters.")
assert_true(!anyNA(analysis_data$assigned_dist_group), "Original distance assignment is missing.")
assert_true(!anyNA(analysis_data$dewormed), "Final take-up is missing.")
assert_true(
  !anyNA(analysis_data %>% select(female, age.census, mu_d, county, assigned_treatment)),
  "A regression covariate is missing."
)
assert_true(
  all(!is.na(analysis_data$dewormed_day[analysis_data$dewormed])),
  "A directly recorded taker has no campaign day."
)
assert_true(
  all(analysis_data$dewormed_day[analysis_data$dewormed] %in% 1:12),
  "A directly recorded take-up day falls outside days 1--12."
)

days <- 1:12
for (day in days) {
  analysis_data[[sprintf("cumulative_%02d", day)]] <-
    as.integer(analysis_data$dewormed & analysis_data$dewormed_day <= day)
  analysis_data[[sprintf("incidence_%02d", day)]] <-
    as.integer(analysis_data$dewormed & analysis_data$dewormed_day == day)
}

cumulative_names <- sprintf("cumulative_%02d", days)
incidence_names <- sprintf("incidence_%02d", days)
outcome_names <- c(cumulative_names, incidence_names)

assert_true(
  all(rowSums(as.matrix(analysis_data[cumulative_names])[, -1, drop = FALSE] <
                as.matrix(analysis_data[cumulative_names])[, -length(days), drop = FALSE]) == 0),
  "Cumulative take-up is not weakly increasing."
)
assert_true(
  identical(analysis_data[["cumulative_12"]], as.integer(analysis_data$dewormed)),
  "Day-12 cumulative take-up does not equal the original endpoint."
)
assert_true(
  identical(rowSums(analysis_data[incidence_names]), as.numeric(analysis_data$dewormed)),
  "Daily incidence does not sum to final take-up."
)

design_formula <- ~ assigned_treatment * assigned_dist_group +
  female + age.census + mu_d + county
X <- model.matrix(design_formula, data = analysis_data)
Y <- as.matrix(analysis_data[outcome_names])

cluster_levels <- sort(unique(analysis_data$original_cluster_id))
cluster_index <- match(analysis_data$original_cluster_id, cluster_levels)
assert_true(length(cluster_levels) == 144, "Cluster indexing failed.")

arm_order <- c("control", "calendar", "ink", "bracelet")
distance_order <- c("close", "far")
targets <- crossing(
  assigned_dist_group = distance_order,
  assigned_treatment = arm_order
) %>%
  mutate(
    assigned_treatment = factor(assigned_treatment, levels = levels(analysis_data$assigned_treatment)),
    assigned_dist_group = factor(assigned_dist_group, levels = levels(analysis_data$assigned_dist_group))
  )

# Average each counterfactual design over the covariate distribution observed in
# its randomized distance cell. This reproduces the APE construction used by
# R/reduced-form/functions.R.
P <- map_dfr(seq_len(nrow(targets)), function(i) {
  distance_i <- as.character(targets$assigned_dist_group[[i]])
  arm_i <- as.character(targets$assigned_treatment[[i]])
  new_data <- analysis_data %>%
    filter(as.character(assigned_dist_group) == distance_i) %>%
    mutate(assigned_treatment = factor(arm_i, levels = levels(analysis_data$assigned_treatment)))
  colMeans(model.matrix(design_formula, data = new_data)) %>%
    as.list() %>%
    as_tibble()
}) %>%
  as.matrix()
colnames(P) <- colnames(X)

fit_all_outcomes <- function(observation_weights) {
  XtWX <- crossprod(X, X * observation_weights)
  XtWY <- crossprod(X, Y * observation_weights)
  beta <- solve(XtWX, XtWY)
  P %*% beta
}

message("Estimating point predictions")
point_predictions <- fit_all_outcomes(rep(1, nrow(analysis_data)))

message("Running ", B, " paired cluster-weight bootstrap draws")
set.seed(bootstrap_seed)
cluster_weight_draws <- matrix(
  rgamma(B * length(cluster_levels), shape = 1, rate = 1),
  nrow = B,
  ncol = length(cluster_levels)
)
draw_array <- array(NA_real_, dim = c(B, nrow(targets), length(outcome_names)))
failed_draws <- integer()

for (b in seq_len(B)) {
  draw_result <- tryCatch(
    fit_all_outcomes(cluster_weight_draws[b, cluster_index]),
    error = function(e) NULL
  )
  if (is.null(draw_result)) {
    failed_draws <- c(failed_draws, b)
  } else {
    draw_array[b, , ] <- draw_result
  }
  if (b %% 100 == 0) message("  completed draw ", b, " of ", B)
}

successful <- which(!apply(is.na(draw_array), 1, any))
assert_true(length(successful) > 0.95 * B, "More than five percent of bootstrap draws failed.")
draw_array <- draw_array[successful, , , drop = FALSE]

outcome_map <- tibble(
  outcome_index = seq_along(outcome_names),
  outcome_name = outcome_names,
  outcome_type = rep(c("Cumulative take-up", "Daily incidence"), each = length(days)),
  day = rep(days, 2)
)

summarize_draws <- function(draws, estimate) {
  tibble(
    estimate = estimate,
    std_error = sd(draws),
    conf_low = unname(quantile(draws, 0.025)),
    conf_high = unname(quantile(draws, 0.975)),
    p_value = min(1, 2 * min(mean(draws <= 0), mean(draws >= 0)))
  )
}

target_index <- function(arm, distance) {
  which(
    as.character(targets$assigned_treatment) == arm &
      as.character(targets$assigned_dist_group) == distance
  )
}

level_summary <- map_dfr(seq_len(nrow(targets)), function(target_i) {
  map_dfr(seq_along(outcome_names), function(outcome_i) {
    summarize_draws(
      draw_array[, target_i, outcome_i],
      point_predictions[target_i, outcome_i]
    ) %>%
      mutate(
        estimand = "level",
        assigned_treatment = as.character(targets$assigned_treatment[[target_i]]),
        assigned_dist_group = as.character(targets$assigned_dist_group[[target_i]]),
        contrast = NA_character_,
        outcome_index = outcome_i
      )
  })
})

far_close_summary <- map_dfr(arm_order, function(arm) {
  far_i <- target_index(arm, "far")
  close_i <- target_index(arm, "close")
  map_dfr(seq_along(outcome_names), function(outcome_i) {
    summarize_draws(
      draw_array[, far_i, outcome_i] - draw_array[, close_i, outcome_i],
      point_predictions[far_i, outcome_i] - point_predictions[close_i, outcome_i]
    ) %>%
      mutate(
        estimand = "far_minus_close",
        assigned_treatment = arm,
        assigned_dist_group = "far - close",
        contrast = paste(arm, "Far - Close"),
        outcome_index = outcome_i
      )
  })
})

treatment_effect_summary <- map_dfr(setdiff(arm_order, "control"), function(arm) {
  map_dfr(distance_order, function(distance) {
    arm_i <- target_index(arm, distance)
    control_i <- target_index("control", distance)
    map_dfr(seq_along(outcome_names), function(outcome_i) {
      summarize_draws(
        draw_array[, arm_i, outcome_i] - draw_array[, control_i, outcome_i],
        point_predictions[arm_i, outcome_i] - point_predictions[control_i, outcome_i]
      ) %>%
        mutate(
          estimand = "treatment_minus_control",
          assigned_treatment = arm,
          assigned_dist_group = distance,
          contrast = paste(arm, "- Control"),
          outcome_index = outcome_i
        )
    })
  })
})

interaction_summary <- map_dfr(setdiff(arm_order, "control"), function(arm) {
  arm_far <- target_index(arm, "far")
  arm_close <- target_index(arm, "close")
  control_far <- target_index("control", "far")
  control_close <- target_index("control", "close")
  map_dfr(seq_along(outcome_names), function(outcome_i) {
    draw_contrast <-
      (draw_array[, arm_far, outcome_i] - draw_array[, control_far, outcome_i]) -
      (draw_array[, arm_close, outcome_i] - draw_array[, control_close, outcome_i])
    point_contrast <-
      (point_predictions[arm_far, outcome_i] - point_predictions[control_far, outcome_i]) -
      (point_predictions[arm_close, outcome_i] - point_predictions[control_close, outcome_i])
    summarize_draws(draw_contrast, point_contrast) %>%
      mutate(
        estimand = "treatment_distance_interaction",
        assigned_treatment = arm,
        assigned_dist_group = "far - close",
        contrast = paste0("(", arm, " - Control) Far - Close"),
        outcome_index = outcome_i
      )
  })
})

signal_summary <- map_dfr(seq_along(outcome_names), function(outcome_i) {
  signal_far_draw <- rowMeans(cbind(
    draw_array[, target_index("ink", "far"), outcome_i],
    draw_array[, target_index("bracelet", "far"), outcome_i]
  ))
  signal_close_draw <- rowMeans(cbind(
    draw_array[, target_index("ink", "close"), outcome_i],
    draw_array[, target_index("bracelet", "close"), outcome_i]
  ))
  no_signal_far_draw <- rowMeans(cbind(
    draw_array[, target_index("control", "far"), outcome_i],
    draw_array[, target_index("calendar", "far"), outcome_i]
  ))
  no_signal_close_draw <- rowMeans(cbind(
    draw_array[, target_index("control", "close"), outcome_i],
    draw_array[, target_index("calendar", "close"), outcome_i]
  ))
  signal_far_point <- mean(point_predictions[c(
    target_index("ink", "far"), target_index("bracelet", "far")), outcome_i])
  signal_close_point <- mean(point_predictions[c(
    target_index("ink", "close"), target_index("bracelet", "close")), outcome_i])
  no_signal_far_point <- mean(point_predictions[c(
    target_index("control", "far"), target_index("calendar", "far")), outcome_i])
  no_signal_close_point <- mean(point_predictions[c(
    target_index("control", "close"), target_index("calendar", "close")), outcome_i])

  bind_rows(
    summarize_draws(
      signal_close_draw - no_signal_close_draw,
      signal_close_point - no_signal_close_point
    ) %>% mutate(assigned_dist_group = "close"),
    summarize_draws(
      signal_far_draw - no_signal_far_draw,
      signal_far_point - no_signal_far_point
    ) %>% mutate(assigned_dist_group = "far"),
    summarize_draws(
      (signal_far_draw - no_signal_far_draw) -
        (signal_close_draw - no_signal_close_draw),
      (signal_far_point - no_signal_far_point) -
        (signal_close_point - no_signal_close_point)
    ) %>% mutate(assigned_dist_group = "far - close")
  ) %>%
    mutate(
      estimand = "signal_minus_no_signal",
      assigned_treatment = "signal",
      contrast = "Any Signal - No Signal",
      outcome_index = outcome_i
    )
})

all_summary <- bind_rows(
  level_summary,
  far_close_summary,
  treatment_effect_summary,
  interaction_summary,
  signal_summary
) %>%
  left_join(outcome_map, by = "outcome_index") %>%
  select(
    outcome_type, day, estimand, assigned_treatment, assigned_dist_group,
    contrast, estimate, std_error, conf_low, conf_high, p_value
  ) %>%
  arrange(outcome_type, day, estimand, assigned_dist_group, assigned_treatment)

write_csv(all_summary, file.path(output_dir, "data", "campaign-day-estimates.csv"))
saveRDS(
  list(
    draws = draw_array,
    targets = targets,
    outcomes = outcome_map,
    bootstrap_seed = bootstrap_seed,
    attempted_draws = B,
    successful_draws = length(successful),
    failed_draws = failed_draws
  ),
  file.path(output_dir, "data", "campaign-day-bootstrap-draws.rds"),
  compress = "xz"
)

existing_endpoint <- read_csv(
  "temp-data/tidy-rf-tes/original-distance-itt-tidy-tes.csv",
  show_col_types = FALSE
) %>%
  filter(
    assigned_treatment %in% c("ink", "calendar", "bracelet"),
    assigned_dist_group %in% c("close", "far", "far - close")
  ) %>%
  select(assigned_treatment, assigned_dist_group, existing_estimate = estimate)

generated_endpoint <- all_summary %>%
  filter(
    outcome_type == "Cumulative take-up",
    day == 12,
    estimand %in% c("treatment_minus_control", "treatment_distance_interaction")
  ) %>%
  select(assigned_treatment, assigned_dist_group, generated_estimate = estimate)

endpoint_check <- existing_endpoint %>%
  left_join(
    generated_endpoint,
    by = c("assigned_treatment", "assigned_dist_group"),
    relationship = "one-to-one"
  ) %>%
  mutate(abs_difference = abs(generated_estimate - existing_estimate))

assert_true(!anyNA(endpoint_check$generated_estimate), "Endpoint validation rows are missing.")
assert_true(
  max(endpoint_check$abs_difference) < 1e-8,
  "Day-12 point estimates do not reproduce the original-assignment ATEs."
)
write_csv(endpoint_check, file.path(output_dir, "data", "endpoint-reproduction-check.csv"))

audit <- tibble(
  check = c(
    "analysis_rows", "analysis_clusters", "direct_recorded_takers",
    "campaign_days", "bootstrap_seed", "bootstrap_attempted",
    "bootstrap_successful", "bootstrap_failed", "max_endpoint_abs_difference"
  ),
  value = c(
    nrow(analysis_data), n_distinct(analysis_data$original_cluster_id),
    sum(analysis_data$dewormed), 12, bootstrap_seed, B,
    length(successful), length(failed_draws), max(endpoint_check$abs_difference)
  )
)
write_csv(audit, file.path(output_dir, "data", "campaign-day-audit.csv"))

treatment_labels <- c(
  control = "Control", calendar = "Calendar", ink = "Ink", bracelet = "Bracelet"
)
treatment_colors <- c(
  control = "#222222", calendar = "#0072B2", ink = "#D55E00", bracelet = "#7B3294"
)
treatment_linetypes <- c(control = "solid", calendar = "dashed", ink = "dotdash", bracelet = "longdash")

plot_data <- all_summary %>%
  filter(estimand == "level") %>%
  mutate(
    assigned_treatment = factor(assigned_treatment, levels = arm_order),
    assigned_dist_group = factor(
      assigned_dist_group,
      levels = distance_order,
      labels = c("Close assignment", "Far assignment")
    )
  )

make_path_plot <- function(type, y_label, file_name) {
  data_i <- plot_data %>% filter(outcome_type == type)
  p <- ggplot(
    data_i,
    aes(
      x = day, y = estimate, color = assigned_treatment,
      fill = assigned_treatment, linetype = assigned_treatment,
      group = assigned_treatment
    )
  ) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.08, color = NA) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.1) +
    facet_wrap(~assigned_dist_group, nrow = 1) +
    scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
    scale_color_manual(values = treatment_colors, labels = treatment_labels) +
    scale_fill_manual(values = treatment_colors, labels = treatment_labels) +
    scale_linetype_manual(values = treatment_linetypes, labels = treatment_labels) +
    labs(x = "Campaign day", y = y_label, color = NULL, fill = NULL, linetype = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
  ggsave(file.path(output_dir, "figures", file_name), p, width = 7.1, height = 3.5)
}

make_path_plot("Cumulative take-up", "Adjusted cumulative take-up", "campaign-day-cumulative.pdf")
make_path_plot("Daily incidence", "Adjusted unconditional daily take-up", "campaign-day-incidence.pdf")

format_cell <- function(estimate, low, high) {
  sprintf("\\makecell[c]{%.3f\\\\{[%.3f, %.3f]}}", estimate, low, high)
}

milestones <- c(1, 3, 6, 12)
table_data <- bind_rows(
  all_summary %>%
    filter(outcome_type == "Cumulative take-up", estimand == "level", day %in% milestones) %>%
    mutate(panel = if_else(assigned_dist_group == "close", "A. Close assignment", "B. Far assignment")),
  all_summary %>%
    filter(outcome_type == "Cumulative take-up", estimand == "far_minus_close", day %in% milestones) %>%
    mutate(panel = "C. Far minus Close")
) %>%
  mutate(
    assigned_treatment = factor(assigned_treatment, levels = arm_order),
    cell = format_cell(estimate, conf_low, conf_high),
    panel = factor(panel, levels = c("A. Close assignment", "B. Far assignment", "C. Far minus Close"))
  ) %>%
  select(panel, assigned_treatment, day, cell) %>%
  pivot_wider(names_from = day, values_from = cell, names_prefix = "day_") %>%
  arrange(panel, assigned_treatment)

table_lines <- c(
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & Day 1 & Day 3 & Day 6 & Day 12 \\\\",
  "\\midrule"
)
for (panel_i in levels(table_data$panel)) {
  table_lines <- c(
    table_lines,
    sprintf(
      "\\multicolumn{5}{l}{\\textit{Panel %s}} \\\\",
      str_replace(panel_i, "\\. ", ": ")
    )
  )
  panel_rows <- table_data %>% filter(panel == panel_i)
  for (row_i in seq_len(nrow(panel_rows))) {
    row <- panel_rows[row_i, ]
    table_lines <- c(
      table_lines,
      paste(
        treatment_labels[[as.character(row$assigned_treatment)]],
        row$day_1, row$day_3, row$day_6, row$day_12,
        sep = " & "
      ) %>% paste0(" \\\\")
    )
  }
  if (panel_i != tail(levels(table_data$panel), 1)) table_lines <- c(table_lines, "\\addlinespace")
}
table_lines <- c(table_lines, "\\bottomrule", "\\end{tabular}")
writeLines(table_lines, file.path(output_dir, "tables", "campaign-day-cumulative-table.tex"))

interaction_table_data <- bind_rows(
  all_summary %>%
    filter(
      outcome_type == "Cumulative take-up",
      estimand == "treatment_distance_interaction",
      day %in% milestones
    ) %>%
    mutate(panel = "A. Arm-specific effects"),
  all_summary %>%
    filter(
      outcome_type == "Cumulative take-up",
      estimand == "signal_minus_no_signal",
      assigned_dist_group == "far - close",
      day %in% milestones
    ) %>%
    mutate(panel = "B. Pooled signal effect")
) %>%
  mutate(
    row_label = recode(
      assigned_treatment,
      calendar = "Calendar $\\times$ Far",
      ink = "Ink $\\times$ Far",
      bracelet = "Bracelet $\\times$ Far",
      signal = "Any Signal $\\times$ Far"
    ),
    row_order = match(assigned_treatment, c("calendar", "ink", "bracelet", "signal")),
    cell = format_cell(estimate, conf_low, conf_high),
    panel = factor(panel, levels = c("A. Arm-specific effects", "B. Pooled signal effect"))
  ) %>%
  select(panel, row_label, row_order, day, cell) %>%
  pivot_wider(names_from = day, values_from = cell, names_prefix = "day_") %>%
  arrange(panel, row_order)

interaction_table_lines <- c(
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & Day 1 & Day 3 & Day 6 & Day 12 \\\\",
  "\\midrule"
)
for (panel_i in levels(interaction_table_data$panel)) {
  interaction_table_lines <- c(
    interaction_table_lines,
    sprintf(
      "\\multicolumn{5}{l}{\\textit{Panel %s}} \\\\",
      str_replace(panel_i, "\\. ", ": ")
    )
  )
  panel_rows <- interaction_table_data %>% filter(panel == panel_i)
  for (row_i in seq_len(nrow(panel_rows))) {
    row <- panel_rows[row_i, ]
    interaction_table_lines <- c(
      interaction_table_lines,
      paste(row$row_label, row$day_1, row$day_3, row$day_6, row$day_12, sep = " & ") %>%
        paste0(" \\\\")
    )
  }
  if (panel_i != tail(levels(interaction_table_data$panel), 1)) {
    interaction_table_lines <- c(interaction_table_lines, "\\addlinespace")
  }
}
interaction_table_lines <- c(interaction_table_lines, "\\bottomrule", "\\end{tabular}")
writeLines(
  interaction_table_lines,
  file.path(output_dir, "tables", "campaign-day-interaction-table.tex")
)

milestone_signal <- all_summary %>%
  filter(
    outcome_type == "Cumulative take-up", estimand == "signal_minus_no_signal",
    assigned_dist_group == "far - close", day %in% milestones
  )

signal_sentence <- paste(
  sprintf(
    "day %d: %.3f [%.3f, %.3f]",
    milestone_signal$day, milestone_signal$estimate,
    milestone_signal$conf_low, milestone_signal$conf_high
  ),
  collapse = "; "
)

final_arm_interactions <- all_summary %>%
  filter(
    outcome_type == "Cumulative take-up", day == 12,
    estimand == "treatment_distance_interaction"
  ) %>%
  mutate(label = treatment_labels[assigned_treatment])

arm_sentence <- paste(
  sprintf(
    "%s: %.3f [%.3f, %.3f]",
    final_arm_interactions$label, final_arm_interactions$estimate,
    final_arm_interactions$conf_low, final_arm_interactions$conf_high
  ),
  collapse = "; "
)

results_lines <- c(
  sprintf(
    "The pooled Any Signal--No Signal Far--Close interaction evolves as follows: %s.",
    signal_sentence
  ),
  sprintf(
    "At day 12, the arm-specific treatment-by-Far interactions are %s.",
    arm_sentence
  ),
  sprintf(
    "The day-12 estimates reproduce the existing original-assignment ATEs to a maximum absolute discrepancy of %.2g.",
    max(endpoint_check$abs_difference)
  ),
  sprintf(
    "Inference uses %d successful paired exponential cluster-weight draws out of %d attempted draws.",
    length(successful), B
  )
)
writeLines(results_lines, file.path(output_dir, "campaign-day-generated-results.tex"))

# Keep the paper-facing standalone appendix artifacts synchronized with the
# local review note. The appendix section itself is maintained separately so
# it can be included with docmute like the other robustness modules.
dir.create(file.path(appendix_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(appendix_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
appendix_artifacts <- c(
  file.path(output_dir, "figures", "campaign-day-cumulative.pdf"),
  file.path(output_dir, "figures", "campaign-day-incidence.pdf"),
  file.path(output_dir, "tables", "campaign-day-cumulative-table.tex"),
  file.path(output_dir, "tables", "campaign-day-interaction-table.tex"),
  file.path(output_dir, "campaign-day-generated-results.tex"),
  "presentations/rf-tables/main-specs/rf_original_distance_itt_tbl.tex"
)
appendix_destinations <- c(
  file.path(appendix_dir, "figures", "campaign-day-cumulative.pdf"),
  file.path(appendix_dir, "figures", "campaign-day-incidence.pdf"),
  file.path(appendix_dir, "tables", "campaign-day-cumulative-table.tex"),
  file.path(appendix_dir, "tables", "campaign-day-interaction-table.tex"),
  file.path(appendix_dir, "tables", "campaign-day-generated-results.tex"),
  file.path(appendix_dir, "tables", "rf-original-distance-itt.tex")
)
copy_ok <- map2_lgl(
  appendix_artifacts,
  appendix_destinations,
  ~file.copy(.x, .y, overwrite = TRUE)
)
assert_true(all(copy_ok), "Failed to synchronize campaign-day appendix artifacts.")

message("Campaign-day analysis complete. Outputs written to ", output_dir)
