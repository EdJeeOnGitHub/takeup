#!/usr/bin/env Rscript

# Population-weighted, coverage-preserving policy allocations and resource-cost
# accounting.  This runner deliberately separates the allocation problem
# (minimize sites) from the ex-post break-even calculation.

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")
source("R/policy/cost-sensitivity.R")

parameter_csv <- policy_option_value(
  args, "--parameter-csv",
  "temp-data/policy-cost-sensitivity/fit105-compact.csv"
)
parameter_type <- policy_option_value(args, "--parameter-type", "raw")
analysis_id <- policy_option_value(args, "--analysis-id", "baseline-posterior")
distance_path <- policy_option_value(
  args, "--distance-data", "optim/data/full-many-pots-experiment.rds"
)
output_path <- policy_option_value(
  args, "--output-path", "ref-reports/policy-cost-sensitivity"
)
work_path <- policy_option_value(
  args, "--work-path", "temp-data/policy-cost-sensitivity"
)
cores <- as.integer(policy_option_value(args, "--cores", "8"))
max_draws <- as.integer(policy_option_value(args, "--max-draws", "0"))
solver <- policy_option_value(args, "--solver", "glpk")
include_legacy <- policy_option_value(args, "--include-legacy", "true") == "true"
pooling_rhos <- as.numeric(strsplit(
  policy_option_value(args, "--pooling-rhos", "0,0.5,1"), ",", fixed = TRUE
)[[1L]])
if (any(!is.finite(pooling_rhos)) || any(pooling_rhos < 0 | pooling_rhos > 1)) {
  stop("--pooling-rhos must be comma-separated values between zero and one.")
}
pooling_rhos <- sort(unique(pooling_rhos))

if (!parameter_type %in% c("raw", "canonical")) {
  stop("--parameter-type must be raw or canonical.", call. = FALSE)
}
if (!file.exists(parameter_csv)) stop("Missing parameter CSV: ", parameter_csv)
if (!file.exists(distance_path)) stop("Missing distance data: ", distance_path)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(work_path, recursive = TRUE, showWarnings = FALSE)

raw_parameters <- read.csv(
  parameter_csv, check.names = FALSE, stringsAsFactors = FALSE
)
if (max_draws > 0L) raw_parameters <- head(raw_parameters, max_draws)
if (nrow(raw_parameters) == 0L) stop("Parameter CSV has no rows.", call. = FALSE)

parameters <- if (parameter_type == "canonical") {
  required <- c(
    "beta_control", "beta_ink", "beta_calendar", "beta_bracelet",
    "dist_beta", "mu_control", "mu_ink", "mu_calendar", "mu_bracelet",
    "mu_dist_control", "mu_dist_ink", "mu_dist_calendar",
    "mu_dist_bracelet", "base_mu_rep", "u_sd", "total_error_sd"
  )
  missing <- setdiff(required, names(raw_parameters))
  if (length(missing)) stop(
    "Canonical parameter file is missing: ", paste(missing, collapse = ", ")
  )
  raw_parameters
} else {
  rows <- lapply(seq_len(nrow(raw_parameters)), function(index) {
    source_csv <- if ("source_csv" %in% names(raw_parameters)) {
      raw_parameters$source_csv[index]
    } else parameter_csv
    canonical_policy_parameters(
      raw_parameters[index, , drop = FALSE], index, source_csv
    )
  })
  do.call(rbind, rows)
}
if (!"draw" %in% names(parameters)) parameters$draw <- seq_len(nrow(parameters))
if (!"replicate" %in% names(parameters)) parameters$replicate <- parameters$draw
if (anyDuplicated(parameters$draw)) stop("Draw identifiers are not unique.")

distance_object <- readRDS(distance_path)
if (!identical(distance_object$candidate_site_mode, "all")) {
  stop("Policy distance object is not the canonical all-candidate object.")
}
if (nrow(distance_object$village_df) != 144L ||
    nrow(distance_object$pot_df) != 1451L ||
    length(unique(distance_object$pot_df$cluster.id)) != 1451L) {
  stop("Canonical policy geography failed its 144-village/1,451-site audit.")
}
villages <- distance_object$village_df
edges <- distance_object$long_distance_mat[, c(
  "index_i", "index_j", "dist", "dist_km"
)]
names(edges) <- c("village_i", "pot_j", "distance", "distance_km")
edges <- edges[edges$distance <= 3500, ]
edges <- edges[order(edges$village_i, edges$pot_j), ]
if (!setequal(unique(edges$village_i), seq_len(nrow(villages)))) {
  stop("The 3.5 km cap leaves at least one village without a candidate site.")
}

# Match each experimental PoT to its unique nearest candidate school. Policy
# assignments elsewhere are the pre-specified set on which consolidation can
# attenuate Bracelet observability toward the Control schedule.
candidate_sites <- distance_object$pot_df
local_pot <- vapply(seq_len(nrow(villages)), function(index) {
  which.min(
    (candidate_sites$lon - villages$pot.lon[index])^2 +
      (candidate_sites$lat - villages$pot.lat[index])^2
  )
}, integer(1))
if (anyDuplicated(local_pot)) {
  stop("Experimental PoT-to-candidate matching is not one-to-one.")
}

census_environment <- new.env(parent = emptyenv())
load("data/takeup_census.RData", envir = census_environment)
census <- census_environment$census.data
population_by_cluster <- aggregate(
  census$num.individuals,
  by = list(cluster.id = census$cluster.id), sum, na.rm = TRUE
)
names(population_by_cluster)[2L] <- "population"
population <- population_by_cluster$population[
  match(villages$cluster.id, population_by_cluster$cluster.id)
]
if (anyNA(population) || any(population <= 0)) {
  stop("Census populations do not cover all policy villages.")
}

solve_one <- function(index) {
  tryCatch({
    parameter <- as.list(parameters[index, , drop = FALSE])
    parameter$model_family <- "gaussian"
    parameter$sd_of_dist <- distance_object$sd_of_dist
    experimental_demand <- predict_policy_draw(
      parameter, villages$dist.to.pot,
      policy_scenarios[policy_scenarios$scenario == "control", , drop = FALSE]
    )$demand
    target_population <- weighted.mean(experimental_demand, population)
    target_legacy <- mean(experimental_demand)

    demand <- lapply(c("control", "bracelet"), function(regime) {
      predict_policy_draw(
        parameter, edges$distance,
        policy_scenarios[policy_scenarios$scenario == regime, , drop = FALSE]
      )$demand
    })
    names(demand) <- c("control", "bracelet")

    specifications <- data.frame(
      estimand = c(
        "population", "population",
        if (include_legacy) c("legacy", "legacy") else character()
      ),
      regime = rep(c("control", "bracelet"), 1L + include_legacy),
      pooling_rho = NA_real_,
      stringsAsFactors = FALSE
    )
    if (length(pooling_rhos)) {
      pooling_specs <- do.call(rbind, lapply(pooling_rhos, function(rho) {
        data.frame(
          estimand = sprintf("pooling-%03d", round(100 * rho)),
          regime = c("control", "bracelet"), pooling_rho = rho,
          stringsAsFactors = FALSE
        )
      }))
      specifications <- rbind(specifications, pooling_specs)
    }
    allocation_cache <- new.env(parent = emptyenv())
    summaries <- lapply(seq_len(nrow(specifications)), function(row) {
      estimand <- specifications$estimand[row]
      regime <- specifications$regime[row]
      population_estimand <- estimand == "population" || startsWith(estimand, "pooling-")
      weights <- if (population_estimand) population else rep(1, length(population))
      target <- if (population_estimand) target_population else target_legacy
      scenario_demand <- demand[[regime]]
      rho <- specifications$pooling_rho[row]
      effective_rho <- if (population_estimand && regime == "bracelet" &&
                           is.finite(rho)) rho else 0
      cache_key <- if (population_estimand) {
        if (regime == "control") "population-control" else
          sprintf("population-bracelet-%03d", round(100 * effective_rho))
      } else paste("legacy", regime, sep = "-")
      if (regime == "bracelet" && is.finite(rho) && rho > 0) {
        distance_sd <- edges$distance / parameter$sd_of_dist
        benefit <- parameter$beta_control - parameter$dist_beta * distance_sd
        mu <- (1 - rho) * policy_mu_rep(distance_sd, parameter, "bracelet") +
          rho * policy_mu_rep(distance_sd, parameter, "control")
        attenuated_demand <- 1 - pnorm(
          solve_policy_fixedpoint(benefit, mu, parameter$u_sd) /
            parameter$total_error_sd
        )
        away_from_experimental_pot <-
          edges$pot_j != local_pot[edges$village_i]
        scenario_demand[away_from_experimental_pot] <-
          attenuated_demand[away_from_experimental_pot]
      }
      if (exists(cache_key, envir = allocation_cache, inherits = FALSE)) {
        value <- get(cache_key, envir = allocation_cache, inherits = FALSE)
      } else {
        fit <- policy_cost_solve(
          edges = edges, demand = scenario_demand, population = weights,
          target_rate = target, site_cost = 1, solver = solver,
          work_path = work_path
        )
        value <- fit$summary
        assigned_away <-
          fit$allocation$pot_j != local_pot[fit$allocation$village_i]
        value$assigned_away_population_share <- sum(
          fit$allocation$population[assigned_away]
        ) / sum(fit$allocation$population)
        value$roundtrip_participant_km <- sum(
          fit$allocation$expected_takers * 2 * fit$allocation$distance_km
        )
        value$taker_mean_distance_km <- sum(
          fit$allocation$expected_takers * fit$allocation$distance_km
        ) / sum(fit$allocation$expected_takers)
        assign(cache_key, value, envir = allocation_cache)
      }
      cbind(
        analysis_id = analysis_id, draw = parameters$draw[index],
        replicate = parameters$replicate[index], estimand = estimand,
        regime = regime, pooling_rho = rho, value, stringsAsFactors = FALSE
      )
    })
    list(status = data.frame(
      analysis_id = analysis_id, draw = parameters$draw[index],
      replicate = parameters$replicate[index], status = "complete",
      error = "", stringsAsFactors = FALSE
    ), results = do.call(rbind, summaries))
  }, error = function(error) list(
    status = data.frame(
      analysis_id = analysis_id, draw = parameters$draw[index],
      replicate = parameters$replicate[index], status = "failed",
      error = conditionMessage(error), stringsAsFactors = FALSE
    ), results = NULL
  ))
}

indices <- seq_len(nrow(parameters))
computed <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(indices, solve_one, mc.cores = cores, mc.preschedule = FALSE)
} else lapply(indices, solve_one)
statuses <- do.call(rbind, lapply(computed, `[[`, "status"))
completed <- vapply(computed, function(value) !is.null(value$results), logical(1))
results <- if (any(completed)) {
  do.call(rbind, lapply(computed[completed], `[[`, "results"))
} else stop("Every policy allocation failed; inspect the status CSV.")

status_path <- file.path(output_path, paste0("policy-run-audit-", analysis_id, ".csv"))
result_path <- file.path(output_path, paste0("policy-allocation-draws-", analysis_id, ".csv"))
write.csv(statuses, status_path, row.names = FALSE)
write.csv(results, result_path, row.names = FALSE)

pooling_results <- results[startsWith(results$estimand, "pooling-"), ]
if (nrow(pooling_results)) {
  pooling_control <- pooling_results[pooling_results$regime == "control", ]
  pooling_bracelet <- pooling_results[pooling_results$regime == "bracelet", ]
  pooling_paired <- merge(
    pooling_control, pooling_bracelet,
    by = c(
      "analysis_id", "draw", "replicate", "estimand", "pooling_rho"
    ), suffixes = c("_control", "_bracelet")
  )
  pooling_paired$sites_saved <-
    pooling_paired$sites_control - pooling_paired$sites_bracelet
  write.csv(
    pooling_paired,
    file.path(output_path, paste0("policy-pooling-draws-", analysis_id, ".csv")),
    row.names = FALSE
  )
}

population_results <- results[results$estimand == "population", ]
control <- population_results[population_results$regime == "control", ]
bracelet <- population_results[population_results$regime == "bracelet", ]
paired <- merge(
  control, bracelet, by = c("analysis_id", "draw", "replicate", "estimand"),
  suffixes = c("_control", "_bracelet")
)
paired$sites_saved <- paired$sites_control - paired$sites_bracelet
paired$signal_cost_difference <- 0.20 * paired$expected_takers_bracelet
paired$roundtrip_participant_km_difference <-
  paired$roundtrip_participant_km_bracelet -
  paired$roundtrip_participant_km_control

travel_values <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25)
fixed_site_values <- c(100, 250, 500, 1000)
cost_rows <- lapply(travel_values, function(travel_value) {
  non_site_cost <- paired$signal_cost_difference +
    travel_value * paired$roundtrip_participant_km_difference
  threshold_raw <- ifelse(
    paired$sites_saved > 0, non_site_cost / paired$sites_saved, NA_real_
  )
  base <- data.frame(
    analysis_id = paired$analysis_id, draw = paired$draw,
    replicate = paired$replicate, travel_cost = travel_value,
    sites_saved = paired$sites_saved,
    signal_cost_difference = paired$signal_cost_difference,
    roundtrip_participant_km_difference =
      paired$roundtrip_participant_km_difference,
    break_even_site_cost_raw = threshold_raw,
    break_even_site_cost = pmax(0, threshold_raw),
    stringsAsFactors = FALSE
  )
  for (fixed_cost in fixed_site_values) {
    base[[paste0("cost_saved_at_", fixed_cost)]] <-
      fixed_cost * paired$sites_saved - non_site_cost
  }
  base
})
costs <- do.call(rbind, cost_rows)
write.csv(
  costs,
  file.path(output_path, paste0("policy-break-even-draws-", analysis_id, ".csv")),
  row.names = FALSE
)

audit <- data.frame(
  analysis_id = analysis_id,
  parameter_type = parameter_type,
  parameter_rows = nrow(parameters),
  successful_rows = sum(statuses$status == "complete"),
  failed_rows = sum(statuses$status == "failed"),
  villages = nrow(villages), candidate_sites = nrow(distance_object$pot_df),
  feasible_links_3500m = nrow(edges), census_adults = sum(population),
  max_experimental_distance_km = max(villages$dist.to.pot) / 1000,
  signal_cost_per_taker = 0.20,
  pooling_rhos = paste(pooling_rhos, collapse = ";"),
  travel_cost_grid = paste(travel_values, collapse = ";"),
  fixed_site_cost_grid = paste(fixed_site_values, collapse = ";"),
  stringsAsFactors = FALSE
)
write.csv(
  audit,
  file.path(output_path, paste0("policy-input-audit-", analysis_id, ".csv")),
  row.names = FALSE
)
message(
  "Completed ", sum(statuses$status == "complete"), "/", nrow(statuses),
  " policy draws for ", analysis_id, "."
)
