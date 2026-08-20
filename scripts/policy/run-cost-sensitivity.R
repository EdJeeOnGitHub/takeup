#!/usr/bin/env Rscript

# Local, review-only policy exercise.  It writes only below
# ref-reports/policy-cost-sensitivity and temp-data/policy-cost-sensitivity.

source("R/policy/bootstrap.R")
source("R/policy/cost-sensitivity.R")

output_path <- "ref-reports/policy-cost-sensitivity"
temporary_path <- "temp-data/policy-cost-sensitivity"
compact_path <- file.path(temporary_path, "fit105-compact.csv")
distance_path <- "optim/data/full-many-pots-experiment.rds"

if (!file.exists(compact_path)) {
  stop(
    "Missing ", compact_path, ". Run scripts/policy/extract-cmdstan-draws.py ",
    "on the four local fit105 chains first.", call. = FALSE
  )
}
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(temporary_path, recursive = TRUE, showWarnings = FALSE)

# Posterior-median plug-in parameters.  Full uncertainty is deliberately left
# to the preferred exponential cluster-weighted modes when Midway returns.
draws <- read.csv(compact_path, check.names = FALSE, stringsAsFactors = FALSE)
parameter_columns <- setdiff(names(draws), c("chain", "iteration", "source_csv"))
median_draw <- as.data.frame(as.list(vapply(
  draws[parameter_columns], median, numeric(1), na.rm = TRUE
)), check.names = FALSE)
parameter <- as.list(canonical_policy_parameters(
  median_draw, replicate = 1L, mode_csv = compact_path
))
parameter$model_family <- "gaussian"

distance_object <- readRDS(distance_path)
parameter$sd_of_dist <- distance_object$sd_of_dist
villages <- distance_object$village_df

# Census adult populations.
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
  stop("Census populations do not cover all policy villages.", call. = FALSE)
}

base_edges <- distance_object$long_distance_mat[, c("index_i", "index_j", "dist", "dist_km")]
names(base_edges) <- c("village_i", "pot_j", "distance", "distance_km")
base_edges <- base_edges[order(base_edges$village_i, base_edges$pot_j), ]

# Each experimental PoT is represented by the nearest candidate school.  The
# match is one-to-one across all 144 communities.  A policy assignment to any
# other candidate is therefore a conservative, pre-specified indicator that
# consolidation changes the community's information environment.
candidate_sites <- distance_object$pot_df
local_pot <- vapply(seq_len(nrow(villages)), function(index) {
  which.min(
    (candidate_sites$lon - villages$pot.lon[index])^2 +
      (candidate_sites$lat - villages$pot.lat[index])^2
  )
}, integer(1))
if (anyDuplicated(local_pot)) stop("Experimental PoT-to-candidate match is not unique.")

demand_for_visibility <- function(edges, visibility) {
  prediction <- predict_policy_draw(
    parameter, edges$distance,
    policy_scenarios[policy_scenarios$scenario == visibility, , drop = FALSE]
  )
  prediction$demand
}

pooled_bracelet_demand <- function(edges, rho) {
  distance_sd <- edges$distance / parameter$sd_of_dist
  benefit <- parameter$beta_control - parameter$dist_beta * distance_sd
  mu_bracelet <- policy_mu_rep(distance_sd, parameter, "bracelet")
  mu_control <- policy_mu_rep(distance_sd, parameter, "control")
  mu <- (1 - rho) * mu_bracelet + rho * mu_control
  cutoff <- solve_policy_fixedpoint(benefit, mu, parameter$u_sd)
  1 - pnorm(cutoff / parameter$total_error_sd)
}

# Experimental-allocation target under Control observability.  Population-
# weighted scenarios preserve the predicted fraction of adults treated;
# unweighted scenarios reproduce the legacy village-average estimand.
experimental_prediction <- predict_policy_draw(
  parameter, villages$dist.to.pot,
  policy_scenarios[policy_scenarios$scenario == "control", , drop = FALSE]
)$demand
target_unweighted <- mean(experimental_prediction)
target_population <- weighted.mean(experimental_prediction, population)

scenario_grid <- data.frame(
  scenario = c(
    "legacy-3500", "population-3500", "population-2500",
    "cost-site-signal", "cost-travel-low", "cost-travel-high",
    "pooling-rho-0", "pooling-rho-50", "pooling-rho-100"
  ),
  label = c(
    "Legacy village-weighted benchmark", "Population weighted, 3.5 km cap",
    "Population weighted, 2.5 km cap", "Site and signal costs",
    "Travel cost: $0.10/round-trip km", "Travel cost: $0.25/round-trip km",
    "Consolidation: no attenuation", "Consolidation: 50% attenuation",
    "Consolidation: Control observability"
  ),
  cap_m = c(3500, 3500, 2500, 3500, 3500, 3500, 3500, 3500, 3500),
  population_weighted = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  site_cost_unit = c(1, 1, 1, 100, 100, 100, 1, 1, 1),
  signal_cost_unit = c(0, 0, 0, 0.20, 0.20, 0.20, 0, 0, 0),
  travel_cost_unit = c(0, 0, 0, 0, 0.10, 0.25, 0, 0, 0),
  pooling_rho = c(NA, NA, NA, NA, NA, NA, 0, 0.5, 1),
  stringsAsFactors = FALSE
)

results <- allocations <- list()
counter <- 0L
for (row in seq_len(nrow(scenario_grid))) {
  specification <- scenario_grid[row, ]
  edges <- base_edges[base_edges$distance <= specification$cap_m, ]
  if (!setequal(unique(edges$village_i), villages$id)) {
    stop("Distance cap leaves a village without a feasible site: ", specification$scenario)
  }
  scenario_population <- if (specification$population_weighted) population else rep(1, length(population))
  target_rate <- if (specification$population_weighted) target_population else target_unweighted
  control_demand <- demand_for_visibility(edges, "control")
  bracelet_demand <- demand_for_visibility(edges, "bracelet")
  for (regime in c("control", "bracelet")) {
    counter <- counter + 1L
    demand <- if (regime == "control") control_demand else bracelet_demand
    attenuated_edges <- rep(FALSE, nrow(edges))
    if (regime == "bracelet" && is.finite(specification$pooling_rho) &&
        specification$pooling_rho > 0) {
      attenuated_edges <- edges$pot_j != local_pot[edges$village_i]
      sensitivity_demand <- pooled_bracelet_demand(
        edges, specification$pooling_rho
      )
      demand[attenuated_edges] <- sensitivity_demand[attenuated_edges]
    }
    fit <- policy_cost_solve(
      edges = edges, demand = demand, population = scenario_population,
      target_rate = target_rate, site_cost = specification$site_cost_unit,
      signal_cost_per_taker = if (regime == "bracelet")
        specification$signal_cost_unit else 0,
      travel_cost_per_roundtrip_km = specification$travel_cost_unit,
      pooled_demand = NULL,
      work_path = temporary_path
    )
    fit$summary$attenuated_population_share <- sum(
      fit$allocation$population[
        fit$allocation$pot_j != local_pot[fit$allocation$village_i]
      ]
    ) / sum(fit$allocation$population)
    summary <- cbind(
      specification, regime = regime,
      analysis_population = sum(scenario_population), fit$summary
    )
    results[[counter]] <- summary
    allocation <- fit$allocation
    allocation$observability_attenuated <- if (
        regime == "bracelet" && is.finite(specification$pooling_rho) &&
        specification$pooling_rho > 0) {
      allocation$pot_j != local_pot[allocation$village_i]
    } else {
      rep(FALSE, nrow(allocation))
    }
    allocation$scenario <- specification$scenario
    allocation$regime <- regime
    allocations[[counter]] <- allocation
    message(
      specification$scenario, " / ", regime, ": ", fit$summary$sites,
      " sites; take-up ", round(fit$summary$population_weighted_takeup, 3),
      "; cost $", round(fit$summary$total_cost)
    )
  }
}

results <- do.call(rbind, results)
allocations <- do.call(rbind, allocations)
write.csv(results, file.path(output_path, "policy-cost-sensitivity-results.csv"), row.names = FALSE)
saveRDS(allocations, file.path(temporary_path, "policy-cost-sensitivity-allocations.rds"), compress = FALSE)

wide <- merge(
  results[results$regime == "control", ],
  results[results$regime == "bracelet", ],
  by = c("scenario", "label"), suffixes = c("_control", "_bracelet")
)
wide$sites_saved <- wide$sites_control - wide$sites_bracelet
wide$total_cost_saved <- wide$total_cost_control - wide$total_cost_bracelet
wide$treated_difference <- wide$expected_takers_bracelet - wide$expected_takers_control
wide$non_site_cost_difference <-
  (wide$signal_cost_bracelet + wide$travel_cost_bracelet) -
  (wide$signal_cost_control + wide$travel_cost_control)
wide$break_even_site_cost <- ifelse(
  wide$sites_saved > 0,
  wide$non_site_cost_difference / wide$sites_saved,
  NA_real_
)
wide <- wide[match(scenario_grid$scenario, wide$scenario), ]
write.csv(wide, file.path(output_path, "policy-cost-sensitivity-contrasts.csv"), row.names = FALSE)

format_number <- function(x, digits = 0) formatC(x, format = "f", digits = digits, big.mark = ",")
format_currency <- function(x) {
  ifelse(
    x < 0,
    paste0("-\\$", format_number(abs(x))),
    paste0("\\$", format_number(x))
  )
}
escape_tex <- function(x) {
  replacements <- c("$" = "\\$", "%" = "\\%", "&" = "\\&",
                    "_" = "\\_", "#" = "\\#")
  for (character in names(replacements)) {
    x <- gsub(character, replacements[[character]], x, fixed = TRUE)
  }
  x
}
table_lines <- c(
  "\\begin{tabular}{lrrrrrrr}",
  "\\toprule",
  "Scenario & Control sites & Bracelet sites & Sites saved & Control cost & Bracelet cost & Cost saved & Break-even site cost \\\\",
  "\\midrule"
)
for (row in seq_len(nrow(wide))) {
  table_lines <- c(table_lines, paste0(
    escape_tex(wide$label[row]), " & ", wide$sites_control[row], " & ",
    wide$sites_bracelet[row], " & ", wide$sites_saved[row], " & ",
    format_currency(wide$total_cost_control[row]), " & ",
    format_currency(wide$total_cost_bracelet[row]), " & ",
    format_currency(wide$total_cost_saved[row]), " & ",
    ifelse(
      startsWith(wide$scenario[row], "cost-") &&
        is.finite(wide$break_even_site_cost[row]),
      format_currency(wide$break_even_site_cost[row]), "--"
    ), " \\\\"
  ))
}
table_lines <- c(table_lines, "\\bottomrule", "\\end{tabular}")
writeLines(table_lines, file.path(output_path, "policy-cost-sensitivity-table.tex"))

pdf(file.path(output_path, "policy-cost-sensitivity-figure.pdf"), width = 8.5, height = 5.2)
old_par <- par(mar = c(8, 4.5, 1, 4.5))
on.exit({par(old_par); dev.off()}, add = TRUE)
positions <- seq_len(nrow(wide))
figure_labels <- c(
  "Legacy", "Population, 3.5 km", "Population, 2.5 km",
  "Site + signal", "Travel $0.10", "Travel $0.25",
  "No attenuation", "50% attenuation", "Control observability"
)
plot(
  positions, wide$sites_saved, type = "b", pch = 19, xaxt = "n",
  xlab = "", ylab = "PoTs saved under Bracelet", col = "#2166AC",
  ylim = range(c(0, wide$sites_saved), finite = TRUE)
)
axis(1, at = positions, labels = figure_labels, las = 2, cex.axis = 0.78)
abline(h = 0, col = "grey70")
par(new = TRUE)
plot(
  positions, wide$total_cost_saved, type = "b", pch = 17, xaxt = "n",
  yaxt = "n", xlab = "", ylab = "", col = "#B2182B",
  ylim = range(c(0, wide$total_cost_saved), finite = TRUE)
)
axis(4, col.axis = "#B2182B")
mtext("Total cost saved ($, illustrative)", side = 4, line = 3, col = "#B2182B")
legend(
  "topright", legend = c("PoTs saved", "Total cost saved"),
  col = c("#2166AC", "#B2182B"), pch = c(19, 17), lty = 1, bty = "n"
)

write.csv(data.frame(
  quantity = c(
    "posterior_draws", "census_adults", "villages", "candidate_sites",
    "target_unweighted", "target_population_weighted",
    "bracelet_unit_cost_usd", "illustrative_site_cost_usd",
    "travel_cost_grid_usd_per_roundtrip_km"
  ),
  value = c(
    nrow(draws), sum(population), nrow(villages), nrow(distance_object$pot_df),
    target_unweighted, target_population, 0.20, 100, "0, 0.10, 0.25"
  )
), file.path(output_path, "policy-cost-sensitivity-inputs.csv"), row.names = FALSE)

message("Wrote review-only policy sensitivity outputs to ", output_path, ".")
