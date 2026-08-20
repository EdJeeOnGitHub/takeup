#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(kableExtra)
  library(readr)
  library(stringr)
})

fit_version <- "104"
model_name <- "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
evaluation_distance_m <- 1500L
input_dir <- "temp-data/bayesian-sensitivity"
output_stem <- file.path(
  input_dir,
  sprintf("structural-object-sensitivity-%dm", evaluation_distance_m)
)
arms <- c("Control", "Ink", "Calendar", "Bracelet")
likelihood_blocks <- c("Takeup", "Beliefs", "WTP")
metric_columns <- paste0("mean_", likelihood_blocks)

derived <- read_csv(
  file.path(
    input_dir,
    sprintf(
      "bayesian-derived-sensitivity-fit%s-%s.csv",
      fit_version, model_name
    )
  ),
  show_col_types = FALSE
)
parameters <- read_csv(
  file.path(
    input_dir,
    sprintf(
      "bayesian-parameter-sensitivity-fit%s-%s.csv",
      fit_version, model_name
    )
  ),
  show_col_types = FALSE
)
visibility <- read_csv(
  file.path(
    input_dir,
    sprintf("visibility-sensitivity-%dm.csv", evaluation_distance_m)
  ),
  show_col_types = FALSE
)

make_arm_rows <- function(data, pattern, object_name) {
  selected <- data |>
    filter(
      str_detect(quantity, fixed(pattern)),
      str_detect(quantity, fixed(as.character(evaluation_distance_m)))
    ) |>
    transmute(
      object = object_name,
      row = str_to_title(str_extract(quantity, "(control|ink|calendar|bracelet)")),
      mean_Takeup = mean_takeup,
      mean_Beliefs = mean_beliefs,
      mean_WTP = mean_wtp
    )
  selected |>
    mutate(row = factor(row, levels = arms)) |>
    arrange(row) |>
    mutate(row = as.character(row))
}

social_multiplier <- make_arm_rows(
  derived, "social_multiplier", "Social multiplier"
)
social_image <- make_arm_rows(
  derived, "social_image_return", "Net social-image return"
)
private_value <- parameters |>
  filter(quantity == "psi") |>
  transmute(
    object = "Private value",
    row = NA_character_,
    mean_Takeup = mean_takeup,
    mean_Beliefs = mean_beliefs,
    mean_WTP = mean_wtp
  )
visibility <- visibility |>
  transmute(
    object,
    row,
    mean_Takeup,
    mean_Beliefs,
    mean_WTP
  )

arm_rows <- bind_rows(social_multiplier, social_image, visibility)
overall_rows <- arm_rows |>
  group_by(object) |>
  summarise(
    row = NA_character_,
    across(all_of(metric_columns), ~ sqrt(mean(.x^2))),
    .groups = "drop"
  )

object_order <- c(
  "Social multiplier",
  "Net social-image return",
  "Private value",
  "Visibility"
)

table_data <- bind_rows(overall_rows, arm_rows, private_value) |>
  mutate(
    object = factor(object, levels = object_order),
    row_order = case_when(
      is.na(row) ~ 0L,
      TRUE ~ match(row, arms)
    )
  ) |>
  arrange(object, row_order) |>
  mutate(
    object = as.character(object),
    across(all_of(metric_columns), abs)
  ) |>
  select(-row_order)

display_data <- table_data |>
  transmute(
    object,
    row,
    Takeup = mean_Takeup,
    Beliefs = mean_Beliefs,
    WTP = mean_WTP
  )

write_csv(display_data, paste0(output_stem, ".csv"))

format_cell <- function(value, is_max) {
  value_text <- formatC(value, format = "f", digits = 2)
  cell_spec(value_text, format = "latex", bold = is_max)
}

object_labels <- c(
  `Social multiplier` =
    "Social multiplier, $\\tilde{SM}(\\textrm{treat}_z,d)$",
  `Net social-image return` =
    paste0(
      "Social-image return, ",
      "$\\mu(\\textrm{treat}_z,d)",
      "\\Delta(w^\\ast(\\textrm{treat}_z,d))$"
    ),
  `Private value` =
    "Private valuation (Calendar relative to Bracelet), $\\psi$",
  Visibility =
    "Observability, $p_{\\textrm{Observed}}(\\textrm{treat}_z,d)$"
)

display_table <- display_data |>
  rowwise() |>
  mutate(
    largest = max(c_across(all_of(likelihood_blocks))),
    label = if (is.na(row)) {
      cell_spec(
        unname(object_labels[object]),
        format = "latex",
        italic = TRUE,
        escape = FALSE
      )
    } else {
      paste0("\\hspace{1em}", row)
    },
    across(
      all_of(likelihood_blocks),
      ~ format_cell(.x, abs(.x - largest) < .Machine$double.eps^0.5)
    )
  ) |>
  ungroup() |>
  select(label, all_of(likelihood_blocks))

custom_save_latex_table <- function(table, table_filename) {
  table_connection <- file(table_filename)
  on.exit(close(table_connection), add = TRUE)

  clean_table <- table |>
    str_remove(fixed("\\begin{table}")) |>
    str_remove("\\\\caption\\{.*\\}") |>
    str_remove("\\\\label\\{.*\\}") |>
    str_remove("\\\\end\\{table\\}")
  clean_table <- clean_table[nzchar(clean_table)]
  writeLines(clean_table, table_connection)
  invisible(table)
}

latex_table <- display_table |>
  kbl(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    linesep = "",
    col.names = c(
      "Structural object",
      "Deworming take-up",
      "Observability",
      "WTP"
    ),
    align = "lrrr"
  ) |>
  column_spec(1, width = "5.4cm") |>
  row_spec(
    c(5, 10, 11),
    extra_latex_after = "\\midrule"
  )

custom_save_latex_table(latex_table, paste0(output_stem, ".tex"))

message("Wrote ", paste0(output_stem, ".csv"))
message("Wrote ", paste0(output_stem, ".tex"))
