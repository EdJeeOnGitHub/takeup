#!/usr/bin/env Rscript

source("R/policy/bootstrap.R")
args <- commandArgs(trailingOnly = TRUE)
input_path <- policy_option_value(
  args, "--input-path", "ref-reports/policy-cost-sensitivity"
)
distance_path <- policy_option_value(
  args, "--distance-data", "optim/data/full-many-pots-experiment.rds"
)

distance_object <- readRDS(distance_path)
stopifnot(
  identical(distance_object$candidate_site_mode, "all"),
  nrow(distance_object$village_df) == 144L,
  nrow(distance_object$pot_df) == 1451L,
  length(unique(distance_object$pot_df$cluster.id)) == 1451L,
  nrow(distance_object$long_distance_mat) == 144L * 1451L
)

geography <- read.csv(file.path(
  input_path, "policy-distance-cap-diagnostics.csv"
))
stopifnot(
  identical(geography$cap_m, c(2500L, 2750L, 3000L, 3250L, 3500L)),
  all(geography$control_sites >= geography$geographic_minimum_sites),
  all(geography$bracelet_sites >= geography$geographic_minimum_sites),
  geography$control_sites[geography$cap_m == 2500L] ==
    geography$geographic_minimum_sites[geography$cap_m == 2500L],
  geography$bracelet_sites[geography$cap_m == 2500L] ==
    geography$geographic_minimum_sites[geography$cap_m == 2500L]
)

analysis_ids <- c("baseline-posterior", "exponential-cluster-weights")
available <- analysis_ids[file.exists(file.path(
  input_path, paste0("policy-allocation-draws-", analysis_ids, ".csv")
))]
if (!length(available)) stop("No allocation outputs to validate.")

for (analysis_id in available) {
  allocation <- read.csv(file.path(
    input_path, paste0("policy-allocation-draws-", analysis_id, ".csv")
  ))
  audit <- read.csv(file.path(
    input_path, paste0("policy-run-audit-", analysis_id, ".csv")
  ))
  costs <- read.csv(file.path(
    input_path, paste0("policy-break-even-draws-", analysis_id, ".csv")
  ))
  stopifnot(
    all(audit$status == "complete"),
    !anyDuplicated(allocation[c("draw", "estimand", "regime")]),
    all(allocation$target_slack_takers >= -1e-4),
    all(is.finite(allocation$population_weighted_takeup)),
    all(costs$sites_saved == costs$sites_saved[match(
      costs$draw, costs$draw[costs$travel_cost == 0]
    )])
  )
  finite <- is.finite(costs$break_even_site_cost_raw)
  reconstructed <- costs$break_even_site_cost_raw[finite] *
    costs$sites_saved[finite] - costs$signal_cost_difference[finite] -
    costs$travel_cost[finite] *
    costs$roundtrip_participant_km_difference[finite]
  stopifnot(max(abs(reconstructed)) < 1e-7)
}

cat(
  "Validated canonical geography and", length(available),
  "policy inference set(s):", paste(available, collapse = ", "), "\n"
)
