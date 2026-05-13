#!/usr/bin/env Rscript

'Generate survey elicitation format checks table.

Usage:
  create-demand-checks-table.R [--input=<path>] [--output-dir=<dir>] [--output-basename=<name>] [--reference=<path>] [--sms-treatment=<value>]

Options:
  --input=<path>             Clean long knowledge-table RDS input [default: data/clean-data/clean-endline-know-table-data-long.rds].
  --output-dir=<dir>         Directory for generated table files [default: tables].
  --output-basename=<name>   Basename for generated table files [default: demand-checks-table].
  --reference=<path>         Optional CSV with comparison values. If omitted, manuscript values are used.
  --sms-treatment=<value>    Restrict to this sms.treatment value; use "all" for no restriction [default: sms.control].
' -> doc

options(warn = 1)

args <- docopt::docopt(doc)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

input_path <- args[["--input"]]
output_dir <- args[["--output-dir"]]
output_basename <- args[["--output-basename"]]
reference_path <- args[["--reference"]]
sms_treatment <- args[["--sms-treatment"]]

required_columns <- c(
  "PARENT_KEY",
  "know.table.type",
  "num.recognized",
  "dewormed",
  "dewormed.know.only",
  "actual.other.dewormed.any.1",
  "actual.other.dewormed.any.2",
  "sms.treatment"
)

if (!file.exists(input_path)) {
  stop("Input file does not exist: ", input_path, call. = FALSE)
}

df <- readRDS(input_path)
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop("Input is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

if (!identical(sms_treatment, "all")) {
  df <- df %>% filter(.data$sms.treatment == sms_treatment)
}

engaged_respondents <- df %>%
  group_by(.data$PARENT_KEY, .data$know.table.type) %>%
  summarise(
    recognized_any = any(.data$num.recognized > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(.data$recognized_any)

analysis_df <- df %>%
  semi_join(engaged_respondents, by = c("PARENT_KEY", "know.table.type"))

respondent_ns <- engaged_respondents %>%
  count(.data$know.table.type, name = "respondents")

get_n <- function(table_type) {
  n <- respondent_ns$respondents[respondent_ns$know.table.type == table_type]
  if (length(n) == 0) return(NA_integer_)
  as.integer(n)
}

clean_binary <- function(x, yes_value) {
  ifelse(is.na(x), NA_integer_, as.integer(x == yes_value))
}

outcome_data <- list(
  recognition = analysis_df %>%
    transmute(
      PARENT_KEY,
      paired_format = as.integer(.data$know.table.type == "table.B"),
      outcome = if_else(.data$know.table.type == "table.B", .data$num.recognized / 2, .data$num.recognized)
    ),
  dont_know = bind_rows(
    analysis_df %>%
      filter(.data$know.table.type == "table.A", .data$num.recognized == 1) %>%
      transmute(PARENT_KEY, paired_format = 0L, outcome = clean_binary(.data$dewormed, "don't know")),
    analysis_df %>%
      filter(.data$know.table.type == "table.B", .data$num.recognized == 1) %>%
      transmute(PARENT_KEY, paired_format = 1L, outcome = clean_binary(.data$dewormed.know.only, "don't know"))
  ),
  false_positive = bind_rows(
    analysis_df %>%
      filter(
        .data$know.table.type == "table.A",
        .data$num.recognized == 1,
        .data$actual.other.dewormed.any.1 == FALSE
      ) %>%
      transmute(PARENT_KEY, paired_format = 0L, outcome = clean_binary(.data$dewormed, "yes")),
    analysis_df %>%
      filter(
        .data$know.table.type == "table.B",
        .data$num.recognized == 1,
        .data$actual.other.dewormed.any.1 == FALSE,
        .data$actual.other.dewormed.any.2 == FALSE
      ) %>%
      transmute(PARENT_KEY, paired_format = 1L, outcome = clean_binary(.data$dewormed.know.only, "yes"))
  ),
  false_negative = bind_rows(
    analysis_df %>%
      filter(
        .data$know.table.type == "table.A",
        .data$num.recognized == 1,
        .data$actual.other.dewormed.any.1 == TRUE
      ) %>%
      transmute(PARENT_KEY, paired_format = 0L, outcome = clean_binary(.data$dewormed, "no")),
    analysis_df %>%
      filter(
        .data$know.table.type == "table.B",
        .data$num.recognized == 1,
        .data$actual.other.dewormed.any.1 == TRUE,
        .data$actual.other.dewormed.any.2 == TRUE
      ) %>%
      transmute(PARENT_KEY, paired_format = 1L, outcome = clean_binary(.data$dewormed.know.only, "no"))
  )
)

clustered_lm <- function(dat) {
  dat <- dat %>% filter(!is.na(.data$outcome), !is.na(.data$paired_format), !is.na(.data$PARENT_KEY))

  if (requireNamespace("fixest", quietly = TRUE)) {
    fit <- fixest::feols(outcome ~ paired_format, data = dat, cluster = ~PARENT_KEY)
    coefs <- fixest::coeftable(fit)
    return(list(diff = unname(coefs["paired_format", "Estimate"]), p_value = unname(coefs["paired_format", "Pr(>|t|)"])))
  }

  fit <- lm(outcome ~ paired_format, data = dat)
  x <- model.matrix(fit)
  u <- residuals(fit)
  clusters <- dat$PARENT_KEY
  bread <- solve(crossprod(x))
  meat <- matrix(0, ncol = ncol(x), nrow = ncol(x))
  for (cluster in unique(clusters)) {
    idx <- clusters == cluster
    xu <- crossprod(x[idx, , drop = FALSE], u[idx])
    meat <- meat + tcrossprod(xu)
  }
  n <- nrow(x)
  k <- ncol(x)
  g <- length(unique(clusters))
  vcov <- (g / (g - 1)) * ((n - 1) / (n - k)) * bread %*% meat %*% bread
  se <- sqrt(diag(vcov))["paired_format"]
  t_stat <- unname(coef(fit)["paired_format"] / se)
  p_value <- 2 * stats::pt(abs(t_stat), df = g - 1, lower.tail = FALSE)
  list(diff = unname(coef(fit)["paired_format"]), p_value = p_value)
}

clustered_feols <- function(formula, dat, coefficient) {
  formula_vars <- all.vars(formula)
  dat <- dat %>%
    filter(!is.na(.data$PARENT_KEY)) %>%
    filter(stats::complete.cases(across(all_of(formula_vars))))

  if (requireNamespace("fixest", quietly = TRUE)) {
    fit <- fixest::feols(formula, data = dat, cluster = ~PARENT_KEY)
    coefs <- fixest::coeftable(fit)
    return(list(
      estimate = unname(coefs[coefficient, "Estimate"]),
      se = unname(coefs[coefficient, "Std. Error"]),
      p_value = unname(coefs[coefficient, "Pr(>|t|)"])
    ))
  }

  fit <- lm(formula, data = dat)
  x <- model.matrix(fit)
  u <- residuals(fit)
  clusters <- dat$PARENT_KEY
  bread <- solve(crossprod(x))
  meat <- matrix(0, ncol = ncol(x), nrow = ncol(x))
  for (cluster in unique(clusters)) {
    idx <- clusters == cluster
    xu <- crossprod(x[idx, , drop = FALSE], u[idx])
    meat <- meat + tcrossprod(xu)
  }
  n <- nrow(x)
  k <- ncol(x)
  g <- length(unique(clusters))
  vcov <- (g / (g - 1)) * ((n - 1) / (n - k)) * bread %*% meat %*% bread
  se <- sqrt(diag(vcov))[coefficient]
  t_stat <- unname(coef(fit)[coefficient] / se)
  p_value <- 2 * stats::pt(abs(t_stat), df = g - 1, lower.tail = FALSE)
  list(estimate = unname(coef(fit)[coefficient]), se = unname(se), p_value = p_value)
}

summarise_outcome <- function(dat) {
  dat <- dat %>% filter(!is.na(.data$outcome))
  estimate <- clustered_lm(dat)
  tibble(
    table_a = mean(dat$outcome[dat$paired_format == 0]),
    table_b = mean(dat$outcome[dat$paired_format == 1]),
    difference = estimate$diff,
    p_value = estimate$p_value
  )
}

recognition_truth_data <- bind_rows(
  analysis_df %>%
    filter(.data$know.table.type == "table.A", !is.na(.data$actual.other.dewormed.any.1)) %>%
    transmute(
      PARENT_KEY,
      paired_format = 0L,
      true_participant = as.integer(.data$actual.other.dewormed.any.1),
      recognition_per_name = as.numeric(.data$num.recognized)
    ),
  analysis_df %>%
    filter(
      .data$know.table.type == "table.B",
      !is.na(.data$actual.other.dewormed.any.1),
      !is.na(.data$actual.other.dewormed.any.2),
      .data$actual.other.dewormed.any.1 == .data$actual.other.dewormed.any.2
    ) %>%
    transmute(
      PARENT_KEY,
      paired_format = 1L,
      true_participant = as.integer(.data$actual.other.dewormed.any.1),
      recognition_per_name = as.numeric(.data$num.recognized) / 2
    )
) %>%
  filter(!is.na(.data$recognition_per_name), !is.na(.data$true_participant))

recognition_interaction <- clustered_feols(
  recognition_per_name ~ paired_format * true_participant,
  recognition_truth_data,
  "paired_format:true_participant"
)

recognition_truth_cells <- recognition_truth_data %>%
  group_by(.data$paired_format, .data$true_participant) %>%
  summarise(
    mean = mean(.data$recognition_per_name),
    prompts = n(),
    respondents = n_distinct(.data$PARENT_KEY),
    .groups = "drop"
  )

get_truth_cell <- function(paired_format, true_participant, column) {
  value <- recognition_truth_cells %>%
    filter(.data$paired_format == !!paired_format, .data$true_participant == !!true_participant) %>%
    pull(.data[[column]])
  if (length(value) == 0) return(NA_real_)
  value
}

truth_a_nonparticipant <- get_truth_cell(0L, 0L, "mean")
truth_a_participant <- get_truth_cell(0L, 1L, "mean")
truth_b_nonparticipant <- get_truth_cell(1L, 0L, "mean")
truth_b_participant <- get_truth_cell(1L, 1L, "mean")

truth_rows <- tibble(
  row_id = c("table_a", "table_b", "b_minus_a"),
  label = c("Table A (individual)", "Table B (pair, same truth)", "B--A"),
  true_nonparticipant = c(
    truth_a_nonparticipant,
    truth_b_nonparticipant,
    truth_b_nonparticipant - truth_a_nonparticipant
  ),
  true_participant = c(
    truth_a_participant,
    truth_b_participant,
    truth_b_participant - truth_a_participant
  ),
  participant_gradient = c(
    truth_a_participant - truth_a_nonparticipant,
    truth_b_participant - truth_b_nonparticipant,
    recognition_interaction$estimate
  ),
  p_value = c(NA_real_, NA_real_, recognition_interaction$p_value),
  se = c(NA_real_, NA_real_, recognition_interaction$se)
)

truth_counts <- recognition_truth_cells %>%
  mutate(
    row_id = paste0(if_else(.data$paired_format == 0L, "table_a", "table_b"), "_", if_else(.data$true_participant == 0L, "nonparticipant", "participant"))
  ) %>%
  select(row_id, prompts, respondents)

result_rows <- bind_rows(
  summarise_outcome(outcome_data$recognition) %>% mutate(row_id = "recognition", label = "Recognition per name"),
  summarise_outcome(outcome_data$dont_know) %>% mutate(row_id = "dont_know", label = "Don't know share (recognized)"),
  summarise_outcome(outcome_data$false_positive) %>% mutate(row_id = "false_positive", label = "False positive: P(Yes | Truth = 0)"),
  summarise_outcome(outcome_data$false_negative) %>% mutate(row_id = "false_negative", label = "False negative: P(No | Truth = 1)")
) %>%
  select(row_id, label, table_a, table_b, difference, p_value)

csv_rows <- bind_rows(
  result_rows,
  tibble(
    row_id = "respondents",
    label = "Respondents",
    table_a = get_n("table.A"),
    table_b = get_n("table.B"),
    difference = NA_real_,
    p_value = NA_real_
  )
)

fmt3 <- function(x) {
  if (length(x) == 0) return(character())
  as.character(ifelse(is.na(x), "", sprintf("%.3f", x)))
}

fmt_p <- function(x) {
  if (length(x) == 0) return(character())
  as.character(ifelse(is.na(x), "", ifelse(x < 0.001, "$<0.001$", sprintf("%.3f", x))))
}

fmt_n <- function(x) {
  if (length(x) == 0) return(character())
  as.character(ifelse(is.na(x), "", format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)))
}

custom_save_latex_table <- function(table, table_name, table_output_path = output_dir) {
  table_conn <- file(file.path(table_output_path, paste0(table_name, ".tex")))
  on.exit(close(table_conn), add = TRUE)

  clean_table <- table
  clean_table[1] <- gsub("removeme12345", " ", clean_table[1], fixed = TRUE)
  clean_table <- gsub("removeme12345", " ", clean_table, fixed = TRUE)

  clean_table <- clean_table[
    !grepl("^\\\\begin\\{table\\}", clean_table) &
      !grepl("^\\\\caption\\{.*\\}$", clean_table) &
      !grepl("^\\\\label\\{.*\\}$", clean_table) &
      !grepl("^\\\\end\\{table\\}$", clean_table)
  ]
  clean_table <- gsub("^\\\\centering\\s*$", "", clean_table)
  clean_table <- clean_table[nzchar(clean_table)]

  writeLines(clean_table, table_conn)
  invisible(table)
}

custom_save_latex <- custom_save_latex_table

latex_labels <- c(
  recognition = "Recognition per name",
  dont_know = "Don't know share (recognized)",
  false_positive = "False positive: $P(\\text{Yes}\\mid\\text{Truth}=0)$",
  false_negative = "False negative: $P(\\text{No}\\mid\\text{Truth}=1)$"
)

latex_table <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Survey elicitation format checks for experimenter demand and social desirability}",
  "\\label{tab:demand_checks}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & Table A & Table B & B--A & $p$-value\\\\",
  " & (individual) & (pair) & & \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(result_rows))) {
  row <- result_rows[i, ]
  latex_table <- c(
    latex_table,
    paste0(
      latex_labels[[row$row_id]],
      " & ", fmt3(row$table_a),
      " & ", fmt3(row$table_b),
      " & ", fmt3(row$difference),
      " & ", fmt_p(row$p_value),
      "\\\\"
    )
  )
}

latex_table <- c(
  latex_table,
  "\\midrule",
  paste0("Respondents & ", fmt_n(get_n("table.A")), " & ", fmt_n(get_n("table.B")), " & & \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

get_truth_count <- function(row_id, column) {
  value <- truth_counts %>%
    filter(.data$row_id == !!row_id) %>%
    pull(.data[[column]])
  if (length(value) == 0) return(NA_integer_)
  as.integer(value)
}

truth_csv_rows <- truth_rows %>%
  mutate(
    prompts_true_nonparticipant = c(
      get_truth_count("table_a_nonparticipant", "prompts"),
      get_truth_count("table_b_nonparticipant", "prompts"),
      NA_integer_
    ),
    prompts_true_participant = c(
      get_truth_count("table_a_participant", "prompts"),
      get_truth_count("table_b_participant", "prompts"),
      NA_integer_
    ),
    respondents_true_nonparticipant = c(
      get_truth_count("table_a_nonparticipant", "respondents"),
      get_truth_count("table_b_nonparticipant", "respondents"),
      NA_integer_
    ),
    respondents_true_participant = c(
      get_truth_count("table_a_participant", "respondents"),
      get_truth_count("table_b_participant", "respondents"),
      NA_integer_
    )
  )

truth_latex_table <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Recognition by administrative take-up status across elicitation formats}",
  "\\label{tab:demand_recognition_gradient}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & True non-participant & True participant & Participant gradient & $p$-value\\\\",
  "\\midrule"
)

for (i in seq_len(nrow(truth_rows))) {
  row <- truth_rows[i, ]
  truth_latex_table <- c(
    truth_latex_table,
    paste0(
      row$label,
      " & ", fmt3(row$true_nonparticipant),
      " & ", fmt3(row$true_participant),
      " & ", fmt3(row$participant_gradient),
      " & ", fmt_p(row$p_value),
      "\\\\"
    )
  )
}

truth_latex_table <- c(
  truth_latex_table,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

builtin_reference <- tibble::tribble(
  ~row_id, ~column, ~reference_display, ~reference_numeric,
  "recognition", "table_a", "0.478", 0.478,
  "recognition", "table_b", "0.593", 0.593,
  "recognition", "difference", "0.114", 0.114,
  "recognition", "p_value", "<0.001", 0.001,
  "dont_know", "table_a", "0.204", 0.204,
  "dont_know", "table_b", "0.197", 0.197,
  "dont_know", "difference", "-0.007", -0.007,
  "dont_know", "p_value", "0.610", 0.610,
  "false_positive", "table_a", "0.459", 0.459,
  "false_positive", "table_b", "0.426", 0.426,
  "false_positive", "difference", "-0.033", -0.033,
  "false_positive", "p_value", "0.140", 0.140,
  "false_negative", "table_a", "0.163", 0.163,
  "false_negative", "table_b", "0.219", 0.219,
  "false_negative", "difference", "0.056", 0.056,
  "false_negative", "p_value", "0.005", 0.005,
  "respondents", "table_a", "1,141", 1141,
  "respondents", "table_b", "1,170", 1170
)

truth_builtin_reference <- tibble::tribble(
  ~row_id, ~column, ~reference_display, ~reference_numeric,
  "table_a", "true_nonparticipant", "0.434", 0.434,
  "table_a", "true_participant", "0.525", 0.525,
  "table_a", "participant_gradient", "0.091", 0.091,
  "table_b", "true_nonparticipant", "0.548", 0.548,
  "table_b", "true_participant", "0.644", 0.644,
  "table_b", "participant_gradient", "0.096", 0.096,
  "b_minus_a", "true_nonparticipant", "0.114", 0.114,
  "b_minus_a", "true_participant", "0.119", 0.119,
  "b_minus_a", "participant_gradient", "0.005", 0.005,
  "b_minus_a", "se", "0.016", 0.016,
  "b_minus_a", "p_value", "0.754", 0.754
)

read_reference <- function(path) {
  if (is.null(path)) return(builtin_reference)
  ref <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("row_id", "column")
  if (!all(required %in% names(ref))) {
    stop("Reference CSV must include row_id and column columns.", call. = FALSE)
  }
  if (!"reference_numeric" %in% names(ref)) {
    value_col <- intersect(c("value", "current_value", "reference_value"), names(ref))[1]
    if (is.na(value_col)) {
      stop("Reference CSV must include reference_numeric, value, current_value, or reference_value.", call. = FALSE)
    }
    ref <- ref %>% mutate(reference_numeric = as.numeric(gsub(",", "", .data[[value_col]])))
  }
  if (!"reference_display" %in% names(ref)) {
    ref <- ref %>% mutate(reference_display = as.character(.data$reference_numeric))
  }
  ref %>% select(row_id, column, reference_display, reference_numeric)
}

generated_long <- bind_rows(
  csv_rows %>% transmute(row_id, label, column = "table_a", generated_numeric = .data$table_a, generated_display = if_else(.data$row_id == "respondents", fmt_n(.data$table_a), fmt3(.data$table_a))),
  csv_rows %>% transmute(row_id, label, column = "table_b", generated_numeric = .data$table_b, generated_display = if_else(.data$row_id == "respondents", fmt_n(.data$table_b), fmt3(.data$table_b))),
  csv_rows %>% filter(.data$row_id != "respondents") %>% transmute(row_id, label, column = "difference", generated_numeric = .data$difference, generated_display = fmt3(.data$difference)),
  csv_rows %>% filter(.data$row_id != "respondents") %>% transmute(row_id, label, column = "p_value", generated_numeric = .data$p_value, generated_display = if_else(.data$p_value < 0.001, "<0.001", fmt3(.data$p_value)))
)

comparison <- generated_long %>%
  left_join(read_reference(reference_path), by = c("row_id", "column")) %>%
  mutate(
    difference_from_reference = .data$generated_numeric - .data$reference_numeric,
    matches_reference_display = is.na(.data$reference_display) | .data$generated_display == .data$reference_display
  ) %>%
  arrange(factor(.data$row_id, levels = c("recognition", "dont_know", "false_positive", "false_negative", "respondents")), .data$column)

truth_generated_long <- bind_rows(
  truth_csv_rows %>% transmute(row_id, label, column = "true_nonparticipant", generated_numeric = .data$true_nonparticipant, generated_display = fmt3(.data$true_nonparticipant)),
  truth_csv_rows %>% transmute(row_id, label, column = "true_participant", generated_numeric = .data$true_participant, generated_display = fmt3(.data$true_participant)),
  truth_csv_rows %>% transmute(row_id, label, column = "participant_gradient", generated_numeric = .data$participant_gradient, generated_display = fmt3(.data$participant_gradient)),
  truth_csv_rows %>% filter(!is.na(.data$se)) %>% transmute(row_id, label, column = "se", generated_numeric = .data$se, generated_display = fmt3(.data$se)),
  truth_csv_rows %>% filter(!is.na(.data$p_value)) %>% transmute(row_id, label, column = "p_value", generated_numeric = .data$p_value, generated_display = fmt3(.data$p_value))
)

truth_comparison <- truth_generated_long %>%
  left_join(truth_builtin_reference, by = c("row_id", "column")) %>%
  mutate(
    difference_from_reference = .data$generated_numeric - .data$reference_numeric,
    matches_reference_display = is.na(.data$reference_display) | .data$generated_display == .data$reference_display
  ) %>%
  arrange(factor(.data$row_id, levels = c("table_a", "table_b", "b_minus_a")), .data$column)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- file.path(output_dir, paste0(output_basename, ".csv"))
tex_path <- file.path(output_dir, paste0(output_basename, ".tex"))
comparison_path <- file.path(output_dir, paste0(output_basename, "-comparison.csv"))
truth_output_basename <- "demand-recognition-gradient-table"
truth_csv_path <- file.path(output_dir, paste0(truth_output_basename, ".csv"))
truth_tex_path <- file.path(output_dir, paste0(truth_output_basename, ".tex"))
truth_comparison_path <- file.path(output_dir, paste0(truth_output_basename, "-comparison.csv"))

readr::write_csv(csv_rows, csv_path)
custom_save_latex(latex_table, output_basename)
readr::write_csv(comparison, comparison_path)
readr::write_csv(truth_csv_rows, truth_csv_path)
custom_save_latex(truth_latex_table, truth_output_basename)
readr::write_csv(truth_comparison, truth_comparison_path)

if (any(!comparison$matches_reference_display, na.rm = TRUE)) {
  stop("Generated table differs from reference after formatting; inspect ", comparison_path, call. = FALSE)
}

if (any(!truth_comparison$matches_reference_display, na.rm = TRUE)) {
  stop("Generated recognition-gradient table differs from reference after formatting; inspect ", truth_comparison_path, call. = FALSE)
}

message("Wrote ", tex_path)
message("Wrote ", csv_path)
message("Wrote ", comparison_path)
message("Wrote ", truth_tex_path)
message("Wrote ", truth_csv_path)
message("Wrote ", truth_comparison_path)

# STATA Version Anne used
# * Restrict to sms.control
# keep if sms_treatment == "sms.control"

# * Keep respondents who recognize at least one person/pair in their assigned module
# bysort PARENT_KEY: egen any_recognized = max(num_recognized > 0)
# keep if any_recognized == 1

# * ----------------------------
# * Table A recognition by truth
# * ----------------------------
# preserve
# keep if know_table_type == "table.A"
# keep if !missing(actual_other_dewormed_any_1)

# gen recog_per_name = num_recognized
# gen truth = actual_other_dewormed_any_1

# mean recog_per_name if truth == 1   // 0.525
# mean recog_per_name if truth == 0   // 0.434
# restore

# * ----------------------------
# * Table B recognition by truth
# * same-status pairs only
# * ----------------------------
# preserve
# keep if know_table_type == "table.B"
# keep if !missing(actual_other_dewormed_any_1, actual_other_dewormed_any_2)
# keep if actual_other_dewormed_any_1 == actual_other_dewormed_any_2

# gen recog_per_name = num_recognized / 2
# gen truth = actual_other_dewormed_any_1

# mean recog_per_name if truth == 1   // 0.644
# mean recog_per_name if truth == 0   // 0.548
# restore
# * Build pooled A + B sample
# gen tableB = (know_table_type == "table.B")

# gen recog_per_name = .
# replace recog_per_name = num_recognized if know_table_type == "table.A"
# replace recog_per_name = num_recognized/2 if know_table_type == "table.B"

# gen truth = .
# replace truth = actual_other_dewormed_any_1 if know_table_type == "table.A"
# replace truth = actual_other_dewormed_any_1 if know_table_type == "table.B" ///
#     & actual_other_dewormed_any_1 == actual_other_dewormed_any_2

# keep if !missing(recog_per_name, truth)

# reg recog_per_name i.tableB##i.truth, cluster(PARENT_KEY)