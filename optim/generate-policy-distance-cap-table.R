#!/usr/bin/env Rscript

# Generate the review-only appendix diagnostic showing why a 2.5 km policy cap
# is nearly a no-consolidation constraint under the experimental geography.

source("optim/policy-bootstrap-functions.R")
source("optim/policy-cost-sensitivity-functions.R")

compact_path <- "temp-data/policy-cost-sensitivity/fit105-compact.csv"
distance_path <- "optim/data/full-many-pots-experiment.rds"
csv_path <- "ref-reports/policy-cost-sensitivity/policy-distance-cap-diagnostics.csv"
table_path <- "appendix/structural-robustness/tables/policy-distance-cap-feasibility.tex"
caps <- c(2500, 2750, 3000, 3250, 3500)

if (!file.exists(compact_path)) {
  stop("Missing compact baseline draws: ", compact_path, call. = FALSE)
}

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
if (!identical(distance_object$candidate_site_mode, "all") ||
    nrow(distance_object$pot_df) != 1451L ||
    length(unique(distance_object$pot_df$cluster.id)) != 1451L) {
  stop("Distance-cap table requires the canonical all-1,451-site object.")
}
parameter$sd_of_dist <- distance_object$sd_of_dist
villages <- distance_object$village_df

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

experimental_demand <- predict_policy_draw(
  parameter, villages$dist.to.pot,
  policy_scenarios[policy_scenarios$scenario == "control", , drop = FALSE]
)$demand
target_rate <- weighted.mean(experimental_demand, population)

all_edges <- distance_object$long_distance_mat[, c(
  "index_i", "index_j", "dist", "dist_km"
)]
names(all_edges) <- c("village_i", "pot_j", "distance", "distance_km")

rows <- lapply(caps, function(cap) {
  edges <- all_edges[all_edges$distance <= cap, ]
  edges <- edges[order(edges$village_i, edges$pot_j), ]
  village_degree <- table(factor(edges$village_i, levels = seq_len(nrow(villages))))
  site_degree <- table(factor(
    edges$pot_j, levels = seq_len(nrow(distance_object$pot_df))
  ))

  geographic_floor <- policy_cost_solve(
    edges = edges, demand = rep(1, nrow(edges)),
    population = rep(1, nrow(villages)), target_rate = 0,
    site_cost = 1, work_path = "temp-data/policy-cost-sensitivity"
  )$summary$sites

  sites <- vapply(c("control", "bracelet"), function(regime) {
    demand <- predict_policy_draw(
      parameter, edges$distance,
      policy_scenarios[policy_scenarios$scenario == regime, , drop = FALSE]
    )$demand
    policy_cost_solve(
      edges = edges, demand = demand, population = population,
      target_rate = target_rate, site_cost = 1,
      work_path = "temp-data/policy-cost-sensitivity"
    )$summary$sites
  }, numeric(1))

  data.frame(
    cap_m = cap,
    feasible_links = nrow(edges),
    mean_options_per_village = mean(village_degree),
    villages_with_at_most_two_options = sum(village_degree <= 2),
    sites_feasible_for_multiple_villages = sum(site_degree >= 2),
    geographic_minimum_sites = geographic_floor,
    control_sites = sites["control"],
    bracelet_sites = sites["bracelet"],
    sites_saved = sites["control"] - sites["bracelet"],
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, rows)

dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
write.csv(results, csv_path, row.names = FALSE)

lines <- c(
  "\\begin{tabular}{rrrrrr}",
  "\\toprule",
  "Max. distance & Feasible links & Shareable sites & \\makecell{Geographic\\\\minimum} & \\makecell{Control\\\\sites} & \\makecell{Bracelet\\\\sites} \\\\",
  "\\midrule"
)
for (index in seq_len(nrow(results))) {
  lines <- c(lines, paste0(
    formatC(results$cap_m[index] / 1000, format = "f", digits = 2), " km & ",
    format(results$feasible_links[index], big.mark = ",", scientific = FALSE), " & ",
    results$sites_feasible_for_multiple_villages[index], " & ",
    results$geographic_minimum_sites[index], " & ",
    results$control_sites[index], " & ",
    results$bracelet_sites[index], " \\\\"
  ))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, table_path)
message("Wrote policy distance-cap diagnostic table and CSV.")
