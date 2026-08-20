#!/usr/bin/env Rscript

# Write a custom_save_latex_table-style fragment: no table float, caption,
# label, or notes. Those belong to the LaTeX file that inputs the fragment.
args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")

parameter_path <- policy_option_value(args, "--parameter-csv", paste0(
  "optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/",
  "agg-full-many-pots-exponential-cluster-weights/",
  "policy-bootstrap-parameters.csv"
))
distance_path <- policy_option_value(
  args, "--distance-data", "optim/data/full-many-pots-experiment.rds"
)
validation_path <- policy_option_value(
  args, "--validation-draws",
  "temp-data/main-core-exponential-cluster-weight-999/estimand-draws.csv"
)
draw_path <- policy_option_value(
  args, "--draw-output",
  "temp-data/main-core-exponential-cluster-weight-999/finite-segment-multiplier-draws.csv"
)
summary_path <- policy_option_value(
  args, "--summary-output",
  "temp-data/main-core-exponential-cluster-weight-999/finite-segment-multiplier-summary.csv"
)
table_path <- policy_option_value(
  args, "--table-path",
  paste0("appendix/structural-robustness/tables/",
         "main-core-exponential-cluster-weight-finite-segments.tex")
)

parameters <- read.csv(parameter_path, stringsAsFactors = FALSE)
parameters$sd_of_dist <- readRDS(distance_path)$sd_of_dist
required <- c(
  "replicate", "beta_control", "beta_ink", "beta_calendar",
  "beta_bracelet", "dist_beta", "mu_control", "mu_ink", "mu_calendar",
  "mu_bracelet", "mu_dist_control", "mu_dist_ink", "mu_dist_calendar",
  "mu_dist_bracelet", "base_mu_rep", "u_sd", "sd_of_dist"
)
missing <- setdiff(required, names(parameters))
if (length(missing)) stop("Missing parameters: ", paste(missing, collapse = ", "))
if (anyDuplicated(parameters$replicate)) stop("Duplicate replicate identifiers.")
if (any(!is.finite(as.matrix(parameters[, setdiff(required, "replicate")]))) ) {
  stop("Non-finite canonical parameters.")
}

# Match core_gq_find_fixedpoint_safe in the compact Stan GQ program.
safe_fixedpoint <- function(benefit, mu_rep, u_sd) {
  lower <- -5.5
  upper <- 5.5
  fl <- policy_fixedpoint_residual(lower, benefit, mu_rep, u_sd)
  fu <- policy_fixedpoint_residual(upper, benefit, mu_rep, u_sd)
  midpoint <- -benefit
  if (fl * fu <= 0) {
    for (step in seq_len(50L)) {
      midpoint <- (lower + upper) / 2
      fm <- policy_fixedpoint_residual(midpoint, benefit, mu_rep, u_sd)
      if (fl * fm <= 0) {
        upper <- midpoint
        fu <- fm
      } else {
        lower <- midpoint
        fl <- fm
      }
    }
  } else {
    for (step in seq_len(100L)) {
      residual <- policy_fixedpoint_residual(midpoint, benefit, mu_rep, u_sd)
      midpoint <- min(5.5, max(-5.5, midpoint - 0.5 * residual))
    }
  }
  midpoint
}

arms <- c("Bracelet", "Calendar", "Ink", "Control")
distances <- c(0, 500, 1500, 2500)
cutoffs <- do.call(rbind, lapply(seq_len(nrow(parameters)), function(i) {
  parameter <- parameters[i, , drop = FALSE]
  do.call(rbind, lapply(arms, function(arm) {
    key <- tolower(arm)
    intercept <- parameter[[paste0("beta_", key)]]
    if (arm != "Control") intercept <- parameter$beta_control + intercept
    distance_sd <- distances / parameter$sd_of_dist
    benefit <- intercept - parameter$dist_beta * distance_sd
    mu_rep <- policy_mu_rep(distance_sd, parameter, key)
    data.frame(
      replicate = parameter$replicate, treatment = arm,
      distance = distances,
      cutoff = mapply(safe_fixedpoint, benefit, mu_rep,
                      MoreArgs = list(u_sd = parameter$u_sd)),
      dist_beta = parameter$dist_beta, sd_of_dist = parameter$sd_of_dist,
      stringsAsFactors = FALSE
    )
  }))
}))

finite_change <- function(data, from, to) {
  change <- data$cutoff[data$distance == to] - data$cutoff[data$distance == from]
  change / (data$dist_beta[1] * ((to - from) / data$sd_of_dist[1]))
}
groups <- split(cutoffs, interaction(cutoffs$replicate, cutoffs$treatment))
arm_draws <- do.call(rbind, lapply(groups, function(data) data.frame(
  replicate = data$replicate[1], row = data$treatment[1],
  column = c("Combined", "Close segment", "Far segment", "Validation"),
  value = c(
    finite_change(data, 0, 2500), finite_change(data, 0, 1500),
    finite_change(data, 1500, 2500), finite_change(data, 500, 2500)
  ), stringsAsFactors = FALSE
)))
row.names(arm_draws) <- NULL

# Rounded canonical parameters should reproduce the saved 500--2500m Stan GQ
# to within 0.01 before we trust the newly computed endpoints.
saved <- read.csv(validation_path, stringsAsFactors = FALSE)
saved <- saved[
  saved$specification == "Exponential cluster weights" &
    saved$estimand == "finite_multiplier" & saved$subgroup == "500--2500m",
  c("replicate", "treatment", "value")
]
check <- arm_draws[arm_draws$column == "Validation", ]
names(check)[names(check) == "row"] <- "treatment"
check <- merge(check, saved, by = c("replicate", "treatment"),
               suffixes = c("_reconstructed", "_stan"))
if (nrow(check) != 4L * nrow(parameters)) stop("Incomplete validation match.")
max_error <- max(abs(check$value_reconstructed - check$value_stan))
if (!is.finite(max_error) || max_error > 0.01) {
  stop(sprintf("Stan GQ validation failed: maximum error %.6f.", max_error))
}
arm_draws <- arm_draws[arm_draws$column != "Validation", ]

wide <- reshape(arm_draws, idvar = c("replicate", "column"),
                timevar = "row", direction = "wide")
comparisons <- rbind(
  transform(wide[, c("replicate", "column")],
    row = "Signal - No Signal",
    value = (wide$value.Ink + wide$value.Bracelet -
             wide$value.Control - wide$value.Calendar) / 2),
  transform(wide[, c("replicate", "column")],
    row = "Bracelet - Calendar",
    value = wide$value.Bracelet - wide$value.Calendar)
)
draws <- rbind(arm_draws, comparisons[, names(arm_draws)])
row_order <- c("Bracelet", "Calendar", "Ink", "Control",
               "Signal - No Signal", "Bracelet - Calendar")
column_order <- c("Combined", "Close segment", "Far segment")

summary <- do.call(rbind, lapply(row_order, function(row) do.call(rbind,
  lapply(column_order, function(column) {
    value <- draws$value[draws$row == row & draws$column == column]
    data.frame(row, column, estimate = median(value),
      conf_low = unname(quantile(value, .025)),
      conf_high = unname(quantile(value, .975)),
      directional_probability = if (row %in% c(
        "Signal - No Signal", "Bracelet - Calendar"
      )) mean(value < 0) else NA_real_,
      n_refits = length(value))
  })
)))
draws$row <- factor(draws$row, levels = row_order)
draws$column <- factor(draws$column, levels = column_order)
draws <- draws[order(draws$replicate, draws$row, draws$column), ]

dir.create(dirname(draw_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
write.csv(draws, draw_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)

cell <- function(row, column) {
  x <- summary[summary$row == row & summary$column == column, ]
  sprintf("\\makecell[c]{%.3f\\\\(%.3f, %.3f)}",
          x$estimate, x$conf_low, x$conf_high)
}
probability_cell <- function(row, column) {
  x <- summary[summary$row == row & summary$column == column, ]
  sprintf("%.3f", x$directional_probability)
}
lines <- c(
  "\\centering",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
  "\\begin{tabular}[t]{lccc}", "\\toprule",
  paste0("\\multicolumn{1}{c}{ } & ",
         "\\multicolumn{3}{c}{Cluster weighted-likelihood bootstrap} \\\\"),
  "\\cmidrule(l{3pt}r{3pt}){2-4}",
  " & Combined & Close segment & Far segment \\\\ ",
  "Distance span & 0--2.5 km & 0--1.5 km & 1.5--2.5 km \\\\ ",
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Panel A: Treatment-specific finite multipliers}}\\\\",
  "\\addlinespace"
)
add_rows <- function(lines, rows) {
  for (i in seq_along(rows)) {
    row <- rows[i]
    lines <- c(lines, paste0(row, " & ",
      paste(vapply(column_order, function(column) cell(row, column), ""),
            collapse = " & "), "\\\\"))
    if (i < length(rows)) lines <- c(lines, "\\addlinespace")
  }
  lines
}
lines <- add_rows(lines, row_order[1:4])
lines <- c(lines, "\\addlinespace",
  "\\multicolumn{4}{l}{\\textit{Panel B: Directional probabilities}}\\\\",
  "\\addlinespace")
probability_rows <- c(
  "Signal - No Signal" =
    "$\\Pr(M_{\\mathrm{No\\ Signal}}>M_{\\mathrm{Signal}})$",
  "Bracelet - Calendar" =
    "$\\Pr(M_{\\mathrm{Calendar}}>M_{\\mathrm{Bracelet}})$"
)
for (row in names(probability_rows)) {
  lines <- c(lines, paste0(
    probability_rows[[row]], " & ",
    paste(vapply(
      column_order, function(column) probability_cell(row, column), ""
    ), collapse = " & "), "\\\\"
  ))
  if (row != tail(names(probability_rows), 1L)) {
    lines <- c(lines, "\\addlinespace")
  }
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}%", "}")
writeLines(lines, table_path)

message(sprintf("Validated against Stan GQ (max absolute error %.6f).", max_error))
message("Wrote draws: ", draw_path)
message("Wrote summary: ", summary_path)
message("Wrote LaTeX fragment: ", table_path)
