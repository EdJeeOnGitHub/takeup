#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")

draws_path <- main_core_option_value(
  args, "--draws",
  "temp-data/main-core-cluster-robustness/estimand-draws.csv"
)
table_path <- main_core_option_value(
  args, "--table-path",
  "presentations/tables/fit105/main-core-exponential-cluster-weight-overall-te-table.tex"
)
summary_path <- main_core_option_value(
  args, "--summary-path",
  "temp-data/main-core-cluster-robustness/exponential-cluster-weight-structural-results.csv"
)
method <- main_core_option_value(args, "--method", "exponential")
specification <- switch(
  method,
  exponential = "Exponential cluster weights",
  multinomial = "Cluster bootstrap",
  stop("--method must be exponential or multinomial.", call. = FALSE)
)

draws <- read.csv(draws_path, stringsAsFactors = FALSE)
draws <- draws[draws$specification == specification, ]
if (nrow(draws) == 0L) {
  stop("No draws found for ", specification, ".", call. = FALSE)
}

id <- c("replicate", "subgroup")
ates <- draws[
  draws$estimand == "incentive_ate" &
    draws$treatment %in% c("Bracelet", "Calendar", "Ink"),
  c("replicate", "subgroup", "treatment", "value")
]

treatment_cells <- ates
names(treatment_cells)[names(treatment_cells) == "treatment"] <- "row"
names(treatment_cells)[names(treatment_cells) == "subgroup"] <- "column"

ate_wide <- reshape(
  ates, idvar = id, timevar = "treatment", direction = "wide"
)
contrast_cells <- rbind(
  data.frame(
    replicate = ate_wide$replicate,
    column = ate_wide$subgroup,
    row = "Signal - No Signal",
    # ATEs are relative to Control. Therefore
    # mean(Ink, Bracelet) - mean(Control, Calendar) equals
    # (ATE_Ink + ATE_Bracelet - ATE_Calendar) / 2.
    value = (
      ate_wide$value.Bracelet + ate_wide$value.Ink -
        ate_wide$value.Calendar
    ) / 2
  ),
  data.frame(
    replicate = ate_wide$replicate,
    column = ate_wide$subgroup,
    row = "Bracelet - Calendar",
    value = ate_wide$value.Bracelet - ate_wide$value.Calendar
  )
)

control_cells <- draws[
  draws$estimand == "takeup_level" & draws$treatment == "Control" &
    draws$subgroup %in% c("Combined", "Close", "Far"),
  c("replicate", "subgroup", "value")
]
names(control_cells)[names(control_cells) == "subgroup"] <- "column"
control_cells$row <- "Control mean"
control_difference <- draws[
  draws$estimand == "far_minus_close" & draws$treatment == "Control",
  c("replicate", "value")
]
control_difference$column <- "Far - Close"
control_difference$row <- "Control mean"
control_cells <- rbind(
  control_cells[, c("replicate", "column", "row", "value")],
  control_difference[, c("replicate", "column", "row", "value")]
)

add_far_close <- function(cells) {
  wide <- reshape(
    cells[cells$column %in% c("Close", "Far"), ],
    idvar = c("replicate", "row"), timevar = "column", direction = "wide"
  )
  data.frame(
    replicate = wide$replicate,
    column = "Far - Close",
    row = wide$row,
    value = wide$value.Far - wide$value.Close
  )
}
treatment_cells <- rbind(treatment_cells, add_far_close(treatment_cells))
contrast_cells <- rbind(contrast_cells, add_far_close(contrast_cells))

cells <- rbind(
  treatment_cells[, c("replicate", "column", "row", "value")],
  control_cells,
  contrast_cells[, c("replicate", "column", "row", "value")]
)
row_order <- c(
  "Bracelet", "Calendar", "Ink", "Control mean",
  "Signal - No Signal", "Bracelet - Calendar"
)
column_order <- c("Combined", "Close", "Far", "Far - Close")

summarize_value <- function(value) {
  c(
    estimate = median(value),
    conf_low = unname(quantile(value, 0.025)),
    conf_high = unname(quantile(value, 0.975))
  )
}
summary <- aggregate(
  cells$value,
  list(row = cells$row, column = cells$column),
  summarize_value
)
values <- if (is.list(summary$x)) do.call(rbind, summary$x) else summary$x
summary$x <- NULL
summary <- cbind(summary, values)
summary$row <- factor(summary$row, levels = row_order)
summary$column <- factor(summary$column, levels = column_order)
summary <- summary[order(summary$row, summary$column), ]

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, summary_path, row.names = FALSE)

cell_tex <- function(row, column) {
  value <- summary[summary$row == row & summary$column == column, ]
  if (nrow(value) != 1L) stop("Missing table cell: ", row, "/", column)
  sprintf(
    "\\makecell[c]{%.3f\\\\(%.3f, %.3f)}",
    value$estimate, value$conf_low, value$conf_high
  )
}

lines <- c(
  "\\centering",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
  "\\begin{tabular}[t]{lcccc}",
  "\\toprule",
  paste0(
    "\\multicolumn{1}{c}{ } & ",
    if (method == "exponential") {
      "\\multicolumn{4}{c}{Cluster weighted-likelihood bootstrap} \\\\"
    } else {
      "\\multicolumn{4}{c}{County-stratified cluster bootstrap} \\\\"
    }
  ),
  "\\cmidrule(l{3pt}r{3pt}){2-5}",
  paste0(
    " & Combined & Close & Far & Far - Close \\\\"
  ),
  "Dependent variable: Take-up & (1) & (2) & (3) & (4)\\\\",
  "\\midrule"
)
for (row in row_order) {
  if (row == "Control mean") lines <- c(lines, "\\midrule")
  lines <- c(lines, paste0(
    row, " & ",
    paste(vapply(column_order, function(column) cell_tex(row, column), ""),
          collapse = " & "),
    "\\\\"
  ))
  if (row != tail(row_order, 1L)) lines <- c(lines, "\\addlinespace")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}%", "}")
writeLines(lines, table_path)
message("Wrote ", method, " cluster-weight structural table: ", table_path)
