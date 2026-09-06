#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
new_path <- policy_option_value(args, "--new-path")
old_path <- policy_option_value(args, "--old-path", "ref-reports/policy-cost-sensitivity")
output <- policy_option_value(args, "--output", file.path(new_path, "weighting-change-comparison.csv"))
analyses <- strsplit(policy_option_value(args, "--analyses", "baseline-posterior,exponential-cluster-weights"), ",", fixed = TRUE)[[1]]
rows <- list()
for (analysis in analyses) {
  for (version in c("old-incorrect-population", "corrected-adult-population", "matched-equal-community")) {
    path <- if (version == "old-incorrect-population") old_path else new_path
    estimand <- if (version == "matched-equal-community") "legacy" else "population"
    x <- read.csv(file.path(path, paste0("policy-allocation-draws-", analysis, ".csv")))
    a <- x[x$estimand == estimand & x$regime == "control", ]
    b <- x[x$estimand == estimand & x$regime == "bracelet", ]; b <- b[match(a$draw, b$draw), ]
    stopifnot(nrow(a) > 0, identical(a$draw, b$draw), !anyDuplicated(a$draw), all(is.finite(a$sites - b$sites)))
    saving <- a$sites - b$sites
    interval <- as.numeric(quantile(saving, c(.025,.5,.975)))
    costs <- if (estimand == "population") read.csv(file.path(path, paste0("policy-break-even-draws-", analysis, ".csv"))) else NULL
    threshold <- function(travel) {
      if (is.null(costs)) return(NA_real_)
      median(costs$break_even_site_cost[costs$travel_cost == travel], na.rm = TRUE)
    }
    rows[[length(rows) + 1L]] <- data.frame(analysis_id = analysis, version = version, draws = nrow(a),
      control_sites_median = median(a$sites), bracelet_sites_median = median(b$sites),
      paired_sites_saved_low = interval[1], paired_sites_saved_median = interval[2], paired_sites_saved_high = interval[3],
      probability_sites_saved = mean(saving > 0), break_even_site_cost_no_travel = threshold(0),
      break_even_site_cost_travel_0_10 = threshold(.1),
      matched_draws_with_corrected = version != "old-incorrect-population",
      interpretation = if (version == "old-incorrect-population") {
        "Historical incorrect-weight result; older fit/draw provenance, so change is not a pure weighting effect."
      } else "Matched canonical draws with the same draw-specific experimental-Control target convention.")
  }
}
write.csv(do.call(rbind, rows), output, row.names = FALSE)
