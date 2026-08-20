#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")

input_root <- policy_option_value(
  args, "--input-root",
  "temp-data/policy-model-robustness"
)
baseline_csv <- policy_option_value(
  args, "--baseline-csv",
  "temp-data/policy-exponential-cluster-weights/policy-cluster-bootstrap-summary.csv"
)
output_csv <- policy_option_value(
  args, "--output-csv", "temp-data/policy-model-robustness-summary.csv"
)
table_path <- policy_option_value(
  args, "--table-path",
  "presentations/tables/fit105/optim-policy-model-robustness.tex"
)
contrast_table_path <- policy_option_value(
  args, "--contrast-table-path",
  "presentations/tables/fit105/optim-policy-model-robustness-contrasts.tex"
)

alternative_order <- c(
  "correct-observability",
  "second-order-observability", "grouped-lambda", "arm-lambda",
  "student-t5", "cluster-shock"
)
model_paths <- file.path(input_root, alternative_order, "policy-model-summary.csv")
if (any(!file.exists(model_paths))) {
  stop("Missing model summaries: ", paste(model_paths[!file.exists(model_paths)], collapse = ", "),
       call. = FALSE)
}
alternative <- do.call(rbind, lapply(model_paths, read.csv, stringsAsFactors = FALSE))

model_order <- alternative_order
summary <- alternative
if (file.exists(baseline_csv)) {
  baseline <- read.csv(baseline_csv, stringsAsFactors = FALSE)
  baseline$model_id <- "baseline-exponential-cluster-weights"
  baseline$model_label <- "Baseline model (cluster weighted-likelihood bootstrap)"
  baseline$scenario_label[baseline$scenario_id == 0L] <- "Experimental allocation"
  baseline$draws <- baseline$replicates
  baseline$target_infeasible_share <- ifelse(
    baseline$scenario == "suppress-reputation", 1, 0
  )
  baseline$equilibrium_undefined_share <- 0
  baseline <- baseline[, names(alternative), drop = FALSE]
  summary <- rbind(baseline, alternative)
  model_order <- c("baseline-exponential-cluster-weights", alternative_order)
}
summary$model_id <- factor(summary$model_id, levels = model_order)
summary <- summary[order(summary$model_id, summary$scenario_id), ]
summary$model_id <- as.character(summary$model_id)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, output_csv, row.names = FALSE)

format_cell <- function(estimate, lower, upper, digits, integer = FALSE) {
  if (integer) {
    sprintf("\\makecell[c]{%.0f\\\\(%.0f, %.0f)}", estimate, lower, upper)
  } else {
    format <- paste0("\\makecell[c]{%.", digits, "f\\\\(%.", digits,
                     "f, %.", digits, "f)}")
    sprintf(format, estimate, lower, upper)
  }
}
row_tex <- function(row) {
  diagnostic <- sprintf(
    "%.2f / %.2f", 100 * row$equilibrium_undefined_share,
    100 * row$target_infeasible_share
  )
  paste0(
    row$scenario_label, " & ",
    format_cell(row$n_pot_estimate, row$n_pot_low, row$n_pot_high, 0, TRUE), " & ",
    format_cell(row$takeup_estimate, row$takeup_low, row$takeup_high, 3), " & ",
    format_cell(row$distance_estimate, row$distance_low, row$distance_high, 2), " & ",
    diagnostic, " \\\\"
  )
}

lines <- c(
  "\\small",
  "\\begin{longtable}{p{0.34\\linewidth}cccc}",
  "\\caption{Optimal policy robustness across structural models}\\label{tab:policy-model-robustness}\\\\",
  "\\toprule",
  "Policy information & \\makecell{Assigned\\\\PoTs} & \\makecell{Mean\\\\take-up} & \\makecell{Mean distance\\\\(km)} & \\makecell{Undefined /\\\\infeasible (\\%)} \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{5}{c}{\\tablename\\ \\thetable{} -- continued} \\\\",
  "\\toprule",
  "Policy information & \\makecell{Assigned\\\\PoTs} & \\makecell{Mean\\\\take-up} & \\makecell{Mean distance\\\\(km)} & \\makecell{Undefined /\\\\infeasible (\\%)} \\\\",
  "\\midrule",
  "\\endhead"
)
for (model in model_order) {
  value <- summary[summary$model_id == model, ]
  lines <- c(
    lines,
    paste0("\\multicolumn{5}{l}{\\textit{", value$model_label[1], "}} \\\\"),
    vapply(seq_len(nrow(value)), function(index) row_tex(value[index, ]), character(1)),
    "\\addlinespace"
  )
}
lines <- c(
  lines,
  "\\bottomrule",
  "\\end{longtable}",
  "\\begin{minipage}{\\linewidth}",
  paste0(
    "\\footnotesize \\textit{Notes:} Each panel applies the same five policy counterfactuals, 3.5 km feasible-site cap, and baseline experimental-allocation welfare target. Parentheses contain 2.5th and 97.5th percentiles. ",
    if (file.exists(baseline_csv))
      "The baseline panel uses exponential cluster-weighted mode refits; alternative-model panels use all retained HMC posterior draws. "
    else
      "Displayed panels use all retained HMC posterior draws. ",
    "The final column reports, respectively, the percent of draws with an undefined counterfactual equilibrium and the percent with an unattainable welfare target. Target-infeasible draws report the allocation that maximizes attainable take-up; undefined-equilibrium draws are retained in the diagnostic denominator but excluded from that scenario's quantiles."
  ),
  "\\end{minipage}"
)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, table_path)

contrast_paths <- file.path(
  input_root, alternative_order, "policy-model-contrast-summary.csv"
)
if (any(!file.exists(contrast_paths))) {
  stop("Missing model contrast summaries.", call. = FALSE)
}
contrast <- do.call(rbind, lapply(contrast_paths, read.csv, stringsAsFactors = FALSE))
baseline_replicates_path <- file.path(
  dirname(baseline_csv), "policy-cluster-bootstrap-replicates.csv"
)
if (file.exists(baseline_csv) && file.exists(baseline_replicates_path)) {
  baseline_replicates <- read.csv(
    baseline_replicates_path, stringsAsFactors = FALSE
  )
  baseline_contrast <- do.call(rbind, lapply(list(
    c("control", "bracelet", "Endogenous social image returns"),
    c("static-control", "static-bracelet", "Social image returns fixed at 0.5km")
  ), function(specification) {
    left <- baseline_replicates[baseline_replicates$scenario == specification[1], ]
    right <- baseline_replicates[baseline_replicates$scenario == specification[2], ]
    right <- right[match(left$draw, right$draw), ]
    saved <- left$n_pot - right$n_pot
    interval <- summarize_policy_values(saved)
    data.frame(
      model_id = "baseline-exponential-cluster-weights",
      model_label = "Baseline model (cluster weighted-likelihood bootstrap)",
      contrast = specification[3], left_scenario = specification[1],
      right_scenario = specification[2], left_n_pot = median(left$n_pot),
      right_n_pot = median(right$n_pot), n_pot_saved = interval["estimate"],
      n_pot_saved_low = interval["conf_low"],
      n_pot_saved_high = interval["conf_high"], defined_draws = length(saved),
      draws = length(saved), equilibrium_undefined_share = 0,
      target_infeasible_share = 0, stringsAsFactors = FALSE
    )
  }))
  contrast <- rbind(baseline_contrast, contrast)
}
contrast$model_id <- factor(contrast$model_id, levels = model_order)
contrast <- contrast[order(contrast$model_id, contrast$contrast), ]
contrast$model_id <- as.character(contrast$model_id)
write.csv(contrast, sub("[.]csv$", "-contrasts.csv", output_csv), row.names = FALSE)

dynamic <- contrast[contrast$contrast == "Endogenous social image returns", ]
static <- contrast[contrast$contrast == "Social image returns fixed at 0.5km", ]
static <- static[match(dynamic$model_id, static$model_id), ]
contrast_lines <- c(
  "\\centering",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Structural model & Control PoTs & Bracelet PoTs & PoTs saved & PoTs saved, static \\\\",
  "\\midrule",
  vapply(seq_len(nrow(dynamic)), function(index) {
    saved <- format_cell(
      dynamic$n_pot_saved[index], dynamic$n_pot_saved_low[index],
      dynamic$n_pot_saved_high[index], 0, TRUE
    )
    static_saved <- format_cell(
      static$n_pot_saved[index], static$n_pot_saved_low[index],
      static$n_pot_saved_high[index], 0, TRUE
    )
    paste0(
      dynamic$model_label[index], " & ", sprintf("%.0f", dynamic$left_n_pot[index]),
      " & ", sprintf("%.0f", dynamic$right_n_pot[index]), " & ", saved,
      " & ", static_saved, " \\\\"
    )
  }, character(1)),
  "\\bottomrule",
  "\\end{tabular}}",
  "\\begin{minipage}{0.98\\linewidth}",
  paste0(
    "\\footnotesize \\textit{Notes:} PoTs saved is the draw-level paired difference between the Control and Bracelet optimal allocations at the common welfare target. Parentheses contain 2.5th and 97.5th percentiles. The final column holds social-image returns fixed at their 0.5 km level. ",
    if (file.exists(baseline_replicates_path))
      "The baseline uses exponential cluster-weighted mode refits; other rows use all retained posterior draws. "
    else
      "Rows use all retained posterior draws. ",
    "The companion full table reports undefined-equilibrium and target-infeasibility shares; target-infeasible draws enter paired contrasts at their best feasible allocation."
  ),
  "\\end{minipage}"
)
writeLines(contrast_lines, contrast_table_path)
message("Wrote policy-model robustness summary and table.")
