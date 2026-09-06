#!/usr/bin/env Rscript

# Render the paper-facing optimal-policy robustness table. Each structural
# specification is a row and each policy counterfactual is a column. Cells
# report assigned PoTs and their 95 percent credible intervals.

args <- commandArgs(trailingOnly = TRUE)

option_value <- function(flag, default) {
  prefix <- paste0(flag, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Option supplied more than once: ", flag, call. = FALSE)
  if (length(hit) == 0L) return(default)
  sub(prefix, "", hit, fixed = TRUE)
}

input_root <- option_value(
  "--input-root",
  "temp-data/policy-model-robustness-complete-20260827-002500"
)
output_path <- option_value(
  "--output",
  "appendix/structural-robustness/tables/optim-policy-model-scenarios.tex"
)

model_catalog <- data.frame(
  model_id = c(
    "benchmark", "private-distance-community-image", "full-information",
    "exclude-dispersed", "cluster-shock", "tight-multinomial",
    "second-order-observability", "grouped-lambda", "arm-lambda",
    "student-t5", "finite-mixture"
  ),
  model_label = c(
    "Benchmark", "Individual travel costs",
    "Individual distance observed by peers",
    "Excluding geographically dispersed communities",
    "Unobserved community heterogeneity", "Correct classification of take-up",
    "Perceived community observability", "By public-signal status",
    "By treatment arm", "Heavy-tailed $v$ distribution",
    "Mixture $v$ distribution"
  ),
  family = c(
    "Benchmark specification", "Distance and information assumptions",
    "Distance and information assumptions",
    "Sample and community heterogeneity",
    "Sample and community heterogeneity", "Measurement of observability",
    "Measurement of observability", "Social-image valuation ($\\lambda$)",
    "Social-image valuation ($\\lambda$)",
    "Distribution of intrinsic motivation ($v$)",
    "Distribution of intrinsic motivation ($v$)"
  ),
  model_order = seq_len(11L),
  stringsAsFactors = FALSE
)

if (!dir.exists(input_root)) stop("Missing input root: ", input_root, call. = FALSE)
summary_paths <- file.path(
  input_root, model_catalog$model_id, "policy-model-summary.csv"
)
if (any(!file.exists(summary_paths))) {
  stop(
    "Missing model summaries: ",
    paste(summary_paths[!file.exists(summary_paths)], collapse = ", "),
    call. = FALSE
  )
}
x <- do.call(rbind, lapply(summary_paths, function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}))
required_columns <- c(
  "model_id", "scenario", "n_pot_estimate", "n_pot_low", "n_pot_high",
  "takeup_estimate", "takeup_low", "takeup_high", "distance_estimate",
  "distance_low", "distance_high", "target_infeasible_share"
)
missing_columns <- setdiff(required_columns, names(x))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

scenario_catalog <- data.frame(
  scenario = c("control", "bracelet", "static-control", "static-bracelet", "suppress-reputation"),
  scenario_label = c(
    "Control", "Bracelet", "Control returns at 0.5 km", "Bracelet social image\\\\returns at 0.5 km",
    "No social-image\\\\returns"
  ),
  scenario_order = seq_len(5L),
  stringsAsFactors = FALSE
)

unknown_models <- setdiff(unique(x$model_id), model_catalog$model_id)
if (length(unknown_models)) {
  stop("Add reader-facing labels for model(s): ",
       paste(unknown_models, collapse = ", "), call. = FALSE)
}

# Replace internal labels from upstream summaries with the reader-facing labels
# maintained here.
x$model_label <- NULL
x$scenario_label <- NULL
x <- merge(x, model_catalog, by = "model_id", all.x = TRUE, sort = FALSE)
x <- x[x$scenario %in% scenario_catalog$scenario, , drop = FALSE]
x <- merge(x, scenario_catalog, by = "scenario", all.x = TRUE, sort = FALSE)

duplicate_cells <- duplicated(x[c("model_id", "scenario")])
if (any(duplicate_cells)) {
  stop("Duplicate model-scenario cells in input.", call. = FALSE)
}

models <- unique(x[c("model_id", "model_label", "family", "model_order")])
models <- models[order(models$model_order), , drop = FALSE]
expected <- expand.grid(
  model_id = models$model_id, scenario = scenario_catalog$scenario,
  stringsAsFactors = FALSE
)
observed_key <- paste(x$model_id, x$scenario, sep = "::")
expected_key <- paste(expected$model_id, expected$scenario, sep = "::")
missing_cells <- expected[!expected_key %in% observed_key, , drop = FALSE]
if (nrow(missing_cells)) {
  stop(
    "Missing model-scenario cells: ",
    paste(paste(missing_cells$model_id, missing_cells$scenario, sep = "/"),
          collapse = ", "),
    call. = FALSE
  )
}

fmt <- function(value, digits) sprintf(paste0("%.", digits, "f"), value)
fmt_interval <- function(low, high, digits) {
  paste0("(", fmt(low, digits), ", ", fmt(high, digits), ")")
}
format_cell <- function(row) {
  dagger <- if (is.finite(row$target_infeasible_share) &&
                row$target_infeasible_share > 0) "\\textsuperscript{$\\dagger$}" else ""
  paste0(
    "\\makecell[c]{",
    fmt(row$n_pot_estimate, 0), dagger, "\\\\",
    fmt_interval(row$n_pot_low, row$n_pot_high, 0),
    "}"
  )
}

lines <- c(
  "\\begingroup",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\renewcommand{\\arraystretch}{0.92}",
  "\\begin{longtable}{p{0.27\\linewidth}ccccc}",
  paste0(
    "\\caption{Robustness of optimal policy across structural models}",
    "\\label{tab:policy-model-robustness-scenarios}\\\\"
  ),
  "\\toprule",
  paste0(
    "Structural specification & ",
    paste(paste0("\\makecell[c]{", scenario_catalog$scenario_label, "}"),
          collapse = " & "),
    " \\\\"
  ),
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{6}{c}{\\tablename\\ \\thetable{} -- continued} \\\\",
  "\\toprule",
  paste0(
    "Structural specification & ",
    paste(paste0("\\makecell[c]{", scenario_catalog$scenario_label, "}"),
          collapse = " & "),
    " \\\\"
  ),
  "\\midrule",
  "\\endhead"
)

families <- unique(models$family)
for (family_index in seq_along(families)) {
  family <- families[family_index]
  family_models <- models[models$family == family, , drop = FALSE]
  lines <- c(lines, paste0("\\multicolumn{6}{l}{\\textit{Panel ",
                         LETTERS[family_index], ": ", family, "}} \\\\*"))
  for (index in seq_len(nrow(family_models))) {
    model <- family_models[index, ]
    cells <- vapply(scenario_catalog$scenario, function(scenario_value) {
      row <- x[x$model_id == model$model_id & x$scenario == scenario_value, , drop = FALSE]
      format_cell(row[1, ])
    }, character(1))
    lines <- c(lines, paste0(model$model_label, " & ", paste(cells, collapse = " & "),
                             " \\\\"))
  }
  if (family_index < length(families)) lines <- c(lines, "\\addlinespace")
}

lines <- c(
  lines,
  "\\bottomrule",
  "\\end{longtable}",
  "\\par\\smallskip",
  "\\begin{minipage}{0.98\\linewidth}",
  paste0(
    "\\scriptsize \\textit{Notes:} Each cell reports the posterior median ",
    "number of assigned points of treatment (PoTs); parentheses contain 95 ",
    "percent credible intervals. Each scenario uses the same feasible-site set, 3.5 km ",
    "distance cap, and model/draw-specific experimental Control target with adult-census weights. Returns at 0.5 km ",
    "holds social-image returns fixed at their value at 0.5 km while varying ",
    "distance. \\textsuperscript{$\\dagger$}For at least some posterior draws, ",
    "the common take-up target is infeasible; those draws report the best ",
    "feasible allocation. Undefined counterfactual equilibria are retained in ",
    "diagnostic counts but excluded from the reported quantiles."
  ),
  "\\end{minipage}",
  "\\endgroup"
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, output_path, useBytes = TRUE)
message("Wrote ", output_path)
