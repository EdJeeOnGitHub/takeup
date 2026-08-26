#!/usr/bin/env Rscript

# Render the two paper-facing 1.25 km social-multiplier robustness tables.
#
# The main input is a draw-level CSV with one row per specification, treatment,
# and posterior draw. Required columns:
#
#   table_id       "main" or "prior"
#   spec_id        stable identifier listed in specification_catalog below
#   treatment      Control, Calendar, Ink, or Bracelet
#   draw            posterior draw identifier
#   distance_m      must equal 1250
#   multiplier      reported social multiplier (one = no feedback)
#
# The tight asymmetric-multinomial reporting model is read from the separate
# server-returned exact-1.25 km draw file. It replaces the earlier structural
# binary correct/incorrect sensitivity row; the full model instead treats
# Yes, No, and Don't know as separate reports conditional on administrative
# take-up.
#
# The output is deliberately a custom_save_latex_table-style fragment: it has
# no table float, caption, label, or notes. The paper controls those around an
# \input{...}. The unlabeled third line in each result cell is Pr(M < 1).

args <- commandArgs(trailingOnly = TRUE)

option_value <- function(name, default = NULL) {
  exact <- which(args == name)
  if (length(exact)) {
    if (exact[1] == length(args)) stop("Missing value after ", name)
    return(args[exact[1] + 1L])
  }
  prefix <- paste0(name, "=")
  inline <- startsWith(args, prefix)
  if (any(inline)) return(substring(args[which(inline)[1]], nchar(prefix) + 1L))
  default
}

input_path <- option_value(
  "--input",
  "temp-data/assigned-distance-comparison/multiplier-draws-1250.csv"
)
asymmetric_input_path <- option_value(
  "--asymmetric-input",
  "temp-data/asymmetric-observability-comparison/multiplier-draws-1250.csv"
)
cluster_shock_input_path <- option_value(
  "--cluster-shock-input",
  paste0(
    "temp-data/main-core-cluster-shock-population-gq-1250/",
    "summary/draws.csv"
  )
)
output_dir <- option_value(
  "--output-dir", "appendix/structural-robustness/tables"
)
summary_path <- option_value(
  "--summary-output",
  "temp-data/main-core-multiplier-1250-summary.csv"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(kableExtra)
  library(stringr)
})

reference_distance_m <- 1250
treatment_order <- c("Control", "Calendar", "Ink", "Bracelet")

specification_catalog <- tibble::tribble(
  ~table_id, ~family, ~spec_id, ~specification,
  "main", "Benchmark specification", "baseline", "Benchmark",
  "main", "Distance and information assumptions",
    "private-distance-community-image", "Individual travel costs",
  "main", "Distance and information assumptions", "full-information",
    "Individual distance observed by peers",
  "main", "Sample and community heterogeneity", "exclude-dispersed",
    "Excluding geographically dispersed communities",
  "main", "Sample and community heterogeneity",
    "cluster-random-shock-population", "Unobserved community heterogeneity",
  "main", "Measurement of observability", "tight-multinomial-reporting",
    "Correct classification of take-up",
  "main", "Measurement of observability", "second-order-observability",
    "Perceived community observability",
  "main", "Structural flexibility", "grouped-lambda",
    "Effects varying by public-signal status",
  "main", "Structural flexibility", "arm-specific-lambda",
    "Effects varying by treatment arm",
  "main", "Structural flexibility", "student-t5",
    "Heavy-tailed preferences",
  "main", "Structural flexibility", "finite-mixture",
    "Flexible mixture distribution",
  "prior", "Baseline", "baseline", "Baseline",
  "prior", "Social-image weight", "image-tight", "Tight",
  "prior", "Social-image weight", "image-diffuse", "Diffuse",
  "prior", "Distance cost", "distance-tight", "Tight",
  "prior", "Distance cost", "distance-diffuse", "Diffuse",
  "prior", "Idiosyncratic heterogeneity", "heterogeneity-tight", "Tight",
  "prior", "Idiosyncratic heterogeneity", "heterogeneity-diffuse", "Diffuse",
  "prior", "Observability schedule", "visibility-tight", "Tight",
  "prior", "Observability schedule", "visibility-diffuse", "Diffuse",
  "prior", "Private utility and WTP", "private-tight", "Tight",
  "prior", "Private utility and WTP", "private-diffuse", "Diffuse",
  "prior", "Joint stress test", "joint-tight", "Tight",
  "prior", "Joint stress test", "joint-diffuse", "Diffuse"
) |>
  mutate(spec_order = row_number(), family_order = match(family, unique(family)))

required_columns <- c(
  "table_id", "spec_id", "treatment", "draw", "distance_m", "multiplier"
)
if (!file.exists(input_path)) {
  stop(
    "Missing draw-level input: ", input_path, "\n",
    "Generate or return the exact 1250m multiplier draws before rendering."
  )
}

draws <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
missing_columns <- setdiff(required_columns, names(draws))
if (length(missing_columns)) {
  stop("Input is missing columns: ", paste(missing_columns, collapse = ", "))
}
draws <- draws |>
  select(all_of(required_columns)) |>
  mutate(
    table_id = as.character(table_id),
    spec_id = as.character(spec_id),
    treatment = as.character(treatment),
    distance_m = as.numeric(distance_m),
    multiplier = as.numeric(multiplier)
  ) |>
  filter(.data$spec_id != "correct-classification")

asymmetric_required_columns <- c(
  "specification", "treatment", "draw", "multiplier"
)
if (!file.exists(asymmetric_input_path)) {
  stop(
    "Missing exact-1.25 km asymmetric-reporting draws: ",
    asymmetric_input_path
  )
}
asymmetric_draws <- read.csv(
  asymmetric_input_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
missing_asymmetric_columns <- setdiff(
  asymmetric_required_columns, names(asymmetric_draws)
)
if (length(missing_asymmetric_columns)) {
  stop(
    "Asymmetric-reporting input is missing columns: ",
    paste(missing_asymmetric_columns, collapse = ", ")
  )
}
asymmetric_draws <- asymmetric_draws |>
  filter(.data$specification == "Tight multinomial") |>
  transmute(
    table_id = "main",
    spec_id = "tight-multinomial-reporting",
    treatment = as.character(.data$treatment),
    draw = .data$draw,
    distance_m = reference_distance_m,
    multiplier = as.numeric(.data$multiplier)
  )
if (!nrow(asymmetric_draws)) {
  stop("No Tight multinomial draws found in ", asymmetric_input_path)
}
draws <- bind_rows(draws, asymmetric_draws)

cluster_shock_required_columns <- c(
  "draw", "distance_m", "estimand", "treatment", "value"
)
if (!file.exists(cluster_shock_input_path)) {
  stop(
    "Missing exact-1.25 km cluster-random-shock draws: ",
    cluster_shock_input_path
  )
}
cluster_shock_draws <- read.csv(
  cluster_shock_input_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
missing_cluster_shock_columns <- setdiff(
  cluster_shock_required_columns, names(cluster_shock_draws)
)
if (length(missing_cluster_shock_columns)) {
  stop(
    "Cluster-random-shock input is missing columns: ",
    paste(missing_cluster_shock_columns, collapse = ", ")
  )
}
cluster_shock_draws <- cluster_shock_draws |>
  filter(.data$estimand == "Population-average multiplier") |>
  transmute(
    table_id = "main",
    spec_id = "cluster-random-shock-population",
    treatment = as.character(.data$treatment),
    draw = .data$draw,
    distance_m = as.numeric(.data$distance_m),
    multiplier = as.numeric(.data$value)
  )
if (!nrow(cluster_shock_draws)) {
  stop(
    "No population-average multiplier draws found in ",
    cluster_shock_input_path
  )
}
draws <- bind_rows(draws, cluster_shock_draws)

if (any(!is.finite(draws$distance_m)) ||
    any(draws$distance_m != reference_distance_m)) {
  stop("Every input row must be evaluated at exactly 1250 metres.")
}
if (any(!is.finite(draws$multiplier))) {
  stop("Multiplier draws must all be finite.")
}
if (any(!draws$treatment %in% treatment_order)) {
  stop(
    "Unexpected treatments: ",
    paste(setdiff(unique(draws$treatment), treatment_order), collapse = ", ")
  )
}

expected_keys <- specification_catalog |>
  select(table_id, spec_id) |>
  distinct()
observed_keys <- draws |>
  select(table_id, spec_id) |>
  distinct()
missing_specs <- anti_join(expected_keys, observed_keys, by = c("table_id", "spec_id"))
extra_specs <- anti_join(observed_keys, expected_keys, by = c("table_id", "spec_id"))
if (nrow(missing_specs) || nrow(extra_specs)) {
  message_text <- c()
  if (nrow(missing_specs)) {
    message_text <- c(
      message_text,
      paste0(
        "Missing specifications: ",
        paste(paste(missing_specs$table_id, missing_specs$spec_id, sep = "/"),
              collapse = ", ")
      )
    )
  }
  if (nrow(extra_specs)) {
    message_text <- c(
      message_text,
      paste0(
        "Unexpected specifications: ",
        paste(paste(extra_specs$table_id, extra_specs$spec_id, sep = "/"),
              collapse = ", ")
      )
    )
  }
  stop(paste(message_text, collapse = "\n"))
}

cell_counts <- draws |>
  count(table_id, spec_id, treatment, name = "n_draws")
expected_cells <- tidyr::crossing(
  expected_keys,
  treatment = treatment_order
)
missing_cells <- anti_join(
  expected_cells, cell_counts,
  by = c("table_id", "spec_id", "treatment")
)
if (nrow(missing_cells)) {
  stop(
    "Missing treatment cells: ",
    paste(
      paste(missing_cells$table_id, missing_cells$spec_id,
            missing_cells$treatment, sep = "/"),
      collapse = ", "
    )
  )
}
unbalanced_cells <- cell_counts |>
  group_by(table_id, spec_id) |>
  filter(n_distinct(n_draws) != 1L) |>
  ungroup()
if (nrow(unbalanced_cells)) {
  stop("Treatment arms must contain the same number of draws within a specification.")
}

summary <- draws |>
  group_by(table_id, spec_id, treatment) |>
  summarise(
    median = median(multiplier),
    lower = unname(quantile(multiplier, 0.025)),
    upper = unname(quantile(multiplier, 0.975)),
    probability_mitigation = mean(multiplier < 1),
    n_draws = n(),
    .groups = "drop"
  ) |>
  left_join(
    specification_catalog,
    by = c("table_id", "spec_id"),
    relationship = "many-to-one"
  ) |>
  mutate(
    treatment = factor(treatment, levels = treatment_order),
    distance_m = reference_distance_m
  ) |>
  arrange(spec_order, treatment)

format_number <- function(x) sprintf("%.2f", x)
format_probability <- function(x, n_draws) {
  output <- sprintf("%.3f", x)
  output[x < 0.001] <- "$<0.001$"
  output[x > 0.999] <- "$>0.999$"
  output
}
format_cell <- function(median, lower, upper, probability, n_draws) {
  paste0(
    "\\makecell[c]{", format_number(median),
    "\\\\{[", format_number(lower), ", ", format_number(upper), "]}",
    "\\\\", format_probability(probability, n_draws), "}"
  )
}

custom_save_latex_table <- function(table, table_filename) {
  table_connection <- file(table_filename)
  on.exit(close(table_connection), add = TRUE)
  clean_table <- table |>
    str_remove(fixed("\\begin{table}")) |>
    str_remove("\\\\caption\\{.*\\}") |>
    str_remove("\\\\label\\{.*\\}") |>
    str_remove(fixed("\\end{table}")) |>
    trimws()
  clean_table <- clean_table[nzchar(trimws(clean_table))]
  writeLines(clean_table, table_connection)
  invisible(table)
}

render_table <- function(selected_table_id, filename, caption) {
  catalog <- specification_catalog |>
    filter(.data$table_id == selected_table_id)
  display <- summary |>
    filter(.data$table_id == selected_table_id) |>
    mutate(
      cell = format_cell(
        median, lower, upper, probability_mitigation, n_draws
      )
    ) |>
    select(spec_id, treatment, cell) |>
    tidyr::pivot_wider(names_from = treatment, values_from = cell) |>
    right_join(catalog, by = "spec_id", relationship = "one-to-one") |>
    arrange(spec_order) |>
    select(Specification = specification, all_of(treatment_order))

  latex <- display |>
    kbl(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      linesep = "",
      align = "lcccc",
      caption = caption,
      col.names = c("Specification", treatment_order)
    ) |>
    add_header_above(
      c(" " = 1, "No public signal" = 2, "Public signal" = 2),
      escape = FALSE
    )

  family_runs <- rle(catalog$family)
  ends <- cumsum(family_runs$lengths)
  starts <- ends - family_runs$lengths + 1L
  for (index in rev(seq_along(family_runs$values))) {
    latex <- pack_rows(
      latex,
      group_label = paste0("\\textit{", family_runs$values[index], "}"),
      start_row = starts[index],
      end_row = ends[index],
      escape = FALSE,
      bold = FALSE
    )
  }

  custom_save_latex_table(latex, file.path(output_dir, filename))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

write.csv(
  summary |>
    select(
      table_id, family, spec_id, specification, treatment, distance_m,
      median, lower, upper, probability_mitigation, n_draws
    ),
  summary_path,
  row.names = FALSE
)

render_table(
  "main",
  "main-core-multiplier-robustness-1250m.tex",
  "Social multiplier robustness at 1.25 km"
)
render_table(
  "prior",
  "main-core-multiplier-prior-sensitivity-1250m.tex",
  "Prior sensitivity of the social multiplier at 1.25 km"
)

message("Wrote ", summary_path)
message(
  "Wrote ",
  file.path(output_dir, "main-core-multiplier-robustness-1250m.tex")
)
message(
  "Wrote ",
  file.path(output_dir, "main-core-multiplier-prior-sensitivity-1250m.tex")
)
