#!/usr/bin/env Rscript

# Assemble local-appendix summaries from draw-level population-weighted policy
# allocations.  The exponential cluster-weighted panel is included whenever
# its draw file is present.

args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")

input_path <- policy_option_value(
  args, "--input-path", "ref-reports/policy-cost-sensitivity"
)
table_path <- policy_option_value(
  args, "--table-path", "appendix/structural-robustness/tables"
)
figure_path <- policy_option_value(
  args, "--figure-path", "appendix/structural-robustness/figures"
)
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

analysis_specs <- data.frame(
  analysis_id = c("baseline-posterior", "exponential-cluster-weights"),
  label = c("Baseline posterior", "Cluster weighted-likelihood bootstrap"),
  interval = c("95% credible interval", "95% percentile interval"),
  stringsAsFactors = FALSE
)
allocation_files <- file.path(
  input_path,
  paste0("policy-allocation-draws-", analysis_specs$analysis_id, ".csv")
)
available <- file.exists(allocation_files)
if (!any(available)) stop("No population-policy draw files are available.")
analysis_specs <- analysis_specs[available, , drop = FALSE]
allocation_files <- allocation_files[available]

quantiles <- function(value) {
  as.numeric(quantile(value[is.finite(value)], c(0.025, 0.5, 0.975), names = FALSE))
}
format_number <- function(value, digits = 0) {
  value <- round(value, digits)
  value[value == 0] <- 0
  formatC(value, format = "f", digits = digits, big.mark = ",")
}
format_cell <- function(value, digits = 0, currency = FALSE) {
  values <- quantiles(value)
  formatted <- format_number(values, digits)
  if (currency) formatted <- paste0("\\$", formatted)
  paste0(
    "\\makecell[c]{", formatted[2L], "\\\\{}[", formatted[1L], ", ",
    formatted[3L], "]}"
  )
}

allocation_rows <- contrast_rows <- cost_rows <- pooling_rows <- list()
for (specification in seq_len(nrow(analysis_specs))) {
  analysis_id <- analysis_specs$analysis_id[specification]
  allocations <- read.csv(allocation_files[specification], stringsAsFactors = FALSE)
  for (estimand in intersect(c("legacy", "population"), unique(allocations$estimand))) {
    subset <- allocations[allocations$estimand == estimand, ]
    control <- subset[subset$regime == "control", ]
    bracelet <- subset[subset$regime == "bracelet", ]
    paired <- merge(control, bracelet, by = c("analysis_id", "draw", "replicate"),
                    suffixes = c("_control", "_bracelet"))
    allocation_rows[[length(allocation_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      analysis_label = analysis_specs$label[specification],
      estimand = estimand,
      estimand_label = if (estimand == "population") {
        "Population weighted, 3.5 km"
      } else "Equal village weights, 3.5 km",
      control_sites = format_cell(control$sites),
      bracelet_sites = format_cell(bracelet$sites),
      sites_saved = format_cell(paired$sites_control - paired$sites_bracelet),
      probability_sites_saved = mean(
        paired$sites_control - paired$sites_bracelet > 0
      ),
      control_mean_distance = format_cell(
        control$population_mean_distance_km, 2
      ),
      bracelet_mean_distance = format_cell(
        bracelet$population_mean_distance_km, 2
      ), stringsAsFactors = FALSE
    )
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id, estimand = estimand,
      draw = paired$draw,
      sites_saved = paired$sites_control - paired$sites_bracelet,
      stringsAsFactors = FALSE
    )
  }

  break_even_file <- file.path(
    input_path, paste0("policy-break-even-draws-", analysis_id, ".csv")
  )
  if (!file.exists(break_even_file)) next
  costs <- read.csv(break_even_file, stringsAsFactors = FALSE)
  for (travel_cost in sort(unique(costs$travel_cost))) {
    subset <- costs[costs$travel_cost == travel_cost, ]
    valid_threshold <- subset$break_even_site_cost[is.finite(
      subset$break_even_site_cost
    )]
    cost_rows[[length(cost_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      analysis_label = analysis_specs$label[specification],
      travel_cost = travel_cost,
      break_even = if (length(valid_threshold)) {
        format_cell(valid_threshold, 0, currency = TRUE)
      } else "--",
      median_break_even = if (length(valid_threshold)) median(valid_threshold) else NA_real_,
      lower_break_even = if (length(valid_threshold)) quantiles(valid_threshold)[1L] else NA_real_,
      upper_break_even = if (length(valid_threshold)) quantiles(valid_threshold)[3L] else NA_real_,
      probability_sites_saved = mean(subset$sites_saved > 0),
      probability_cost_saving_100 = mean(subset$cost_saved_at_100 > 0),
      probability_cost_saving_250 = mean(subset$cost_saved_at_250 > 0),
      probability_cost_saving_500 = mean(subset$cost_saved_at_500 > 0),
      probability_cost_saving_1000 = mean(subset$cost_saved_at_1000 > 0),
      stringsAsFactors = FALSE
    )
  }

  pooling_file <- file.path(
    input_path, paste0("policy-pooling-draws-", analysis_id, ".csv")
  )
  if (file.exists(pooling_file)) {
    pooling <- read.csv(pooling_file, stringsAsFactors = FALSE)
    for (rho in sort(unique(pooling$pooling_rho))) {
      subset <- pooling[pooling$pooling_rho == rho, ]
      pooling_rows[[length(pooling_rows) + 1L]] <- data.frame(
        analysis_id = analysis_id,
        analysis_label = analysis_specs$label[specification],
        pooling_rho = rho,
        attenuation_label = if (rho == 0) {
          "No attenuation"
        } else if (rho == 1) {
          "Control observability"
        } else {
          paste0(format_number(100 * rho), "\\% attenuation")
        },
        control_sites = format_cell(subset$sites_control),
        bracelet_sites = format_cell(subset$sites_bracelet),
        sites_saved = format_cell(subset$sites_saved),
        median_sites_saved = median(subset$sites_saved),
        lower_sites_saved = quantiles(subset$sites_saved)[1L],
        upper_sites_saved = quantiles(subset$sites_saved)[3L],
        probability_sites_saved = mean(subset$sites_saved > 0),
        draws = nrow(subset), stringsAsFactors = FALSE
      )
    }
  }
}
allocations <- do.call(rbind, allocation_rows)
costs <- do.call(rbind, cost_rows)
pooling <- if (length(pooling_rows)) do.call(rbind, pooling_rows) else NULL
write.csv(
  allocations,
  file.path(input_path, "policy-population-allocation-summary.csv"),
  row.names = FALSE
)
write.csv(
  costs, file.path(input_path, "policy-break-even-summary.csv"),
  row.names = FALSE
)
if (!is.null(pooling)) {
  write.csv(
    pooling, file.path(input_path, "policy-pooling-summary.csv"),
    row.names = FALSE
  )
}

allocation_tex <- c(
  "\\begin{tabular}{llrrrrrr}", "\\toprule",
  "Inference & Weighting & Control sites & Bracelet sites & Sites saved & $\\Pr(\\text{saved}>0)$ & Control km & Bracelet km \\\\",
  "\\midrule"
)
for (analysis_id in unique(allocations$analysis_id)) {
  block <- allocations[allocations$analysis_id == analysis_id, ]
  for (row in seq_len(nrow(block))) {
    allocation_tex <- c(allocation_tex, paste0(
      if (row == 1L) block$analysis_label[row] else "", " & ",
      block$estimand_label[row], " & ", block$control_sites[row], " & ",
      block$bracelet_sites[row], " & ", block$sites_saved[row], " & ",
      format_number(block$probability_sites_saved[row], 2), " & ",
      block$control_mean_distance[row], " & ",
      block$bracelet_mean_distance[row], " \\\\"
    ))
  }
  if (analysis_id != tail(unique(allocations$analysis_id), 1L)) {
    allocation_tex <- c(allocation_tex, "\\addlinespace")
  }
}
allocation_tex <- c(allocation_tex, "\\bottomrule", "\\end{tabular}")
writeLines(
  allocation_tex,
  file.path(table_path, "policy-population-weighted-allocations.tex")
)

cost_tex <- c(
  "\\begin{tabular}{lrrrrrrr}", "\\toprule",
  "Inference & Travel value & Break-even site cost & $\\Pr(S_B<S_C)$ & $\\Pr(\\Delta C>0)$ at \\$100 & at \\$250 & at \\$500 & at \\$1,000 \\\\",
  "\\midrule"
)
for (analysis_id in unique(costs$analysis_id)) {
  block <- costs[costs$analysis_id == analysis_id, ]
  for (row in seq_len(nrow(block))) {
    cost_tex <- c(cost_tex, paste0(
      if (row == 1L) block$analysis_label[row] else "", " & ",
      "\\$", format_number(block$travel_cost[row], 2), " & ",
      block$break_even[row], " & ",
      format_number(block$probability_sites_saved[row], 2), " & ",
      format_number(block$probability_cost_saving_100[row], 2), " & ",
      format_number(block$probability_cost_saving_250[row], 2), " & ",
      format_number(block$probability_cost_saving_500[row], 2), " & ",
      format_number(block$probability_cost_saving_1000[row], 2), " \\\\"
    ))
  }
  if (analysis_id != tail(unique(costs$analysis_id), 1L)) {
    cost_tex <- c(cost_tex, "\\addlinespace")
  }
}
cost_tex <- c(cost_tex, "\\bottomrule", "\\end{tabular}")
writeLines(cost_tex, file.path(table_path, "policy-break-even-costs.tex"))

if (!is.null(pooling)) {
  pooling_tex <- c(
    "\\begin{tabular}{llrrrr}", "\\toprule",
    "Inference & Consolidation observability & Control sites & Bracelet sites & Sites saved & $\\Pr(\\text{saved}>0)$ \\\\",
    "\\midrule"
  )
  for (analysis_id in unique(pooling$analysis_id)) {
    block <- pooling[pooling$analysis_id == analysis_id, ]
    for (row in seq_len(nrow(block))) {
      pooling_tex <- c(pooling_tex, paste0(
        if (row == 1L) block$analysis_label[row] else "", " & ",
        block$attenuation_label[row], " & ", block$control_sites[row], " & ",
        block$bracelet_sites[row], " & ", block$sites_saved[row], " & ",
        format_number(block$probability_sites_saved[row], 2), " \\\\"
      ))
    }
    if (analysis_id != tail(unique(pooling$analysis_id), 1L)) {
      pooling_tex <- c(pooling_tex, "\\addlinespace")
    }
  }
  pooling_tex <- c(pooling_tex, "\\bottomrule", "\\end{tabular}")
  writeLines(
    pooling_tex,
    file.path(table_path, "policy-pooling-observability.tex")
  )
}

pdf(
  file.path(figure_path, "policy-break-even-site-cost.pdf"),
  width = 7.2, height = 4.8
)
colors <- c("#2166AC", "#B2182B")
plot(
  range(costs$travel_cost), range(c(costs$lower_break_even, costs$upper_break_even),
                                  finite = TRUE),
  type = "n", xlab = "Household travel value ($ per round-trip participant-km)",
  ylab = "Break-even fixed cost per PoT ($)"
)
for (index in seq_along(unique(costs$analysis_id))) {
  analysis_id <- unique(costs$analysis_id)[index]
  block <- costs[costs$analysis_id == analysis_id, ]
  polygon(
    c(block$travel_cost, rev(block$travel_cost)),
    c(block$lower_break_even, rev(block$upper_break_even)),
    col = adjustcolor(colors[index], alpha.f = 0.16), border = NA
  )
  lines(block$travel_cost, block$median_break_even, col = colors[index], lwd = 2)
  points(block$travel_cost, block$median_break_even, col = colors[index], pch = 19)
}
legend(
  "topleft", legend = unique(costs$analysis_label),
  col = colors[seq_along(unique(costs$analysis_id))], lwd = 2, pch = 19,
  bty = "n"
)
dev.off()
message("Wrote population-allocation and break-even appendix artifacts.")
