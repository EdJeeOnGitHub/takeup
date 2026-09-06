#!/usr/bin/env Rscript

# Generate the review-only appendix diagnostic showing why a 2.5 km policy cap
# is nearly a no-consolidation constraint under the experimental geography.

source("R/policy/bootstrap.R")
source("R/policy/cost-sensitivity.R")

args <- commandArgs(trailingOnly = TRUE)
compact_path <- policy_option_value(
  args, "--parameter-csv", "temp-data/policy-cost-sensitivity/fit105-compact.csv"
)
parameter_type <- policy_option_value(args, "--parameter-type", "raw")
distance_path <- policy_option_value(
  args, "--distance-data", "optim/data/full-many-pots-experiment.rds"
)
csv_path <- policy_option_value(
  args, "--csv-path", "ref-reports/policy-cost-sensitivity/policy-distance-cap-diagnostics.csv"
)
table_path <- policy_option_value(
  args, "--table-path", "appendix/structural-robustness/tables/policy-distance-cap-feasibility.tex"
)
caps <- c(2500, 2750, 3000, 3250, 3500)
allocation_path <- policy_option_value(args, "--allocation-path", file.path(dirname(csv_path), "median-allocations"))
canonical_policy_path <- policy_option_value(args, "--policy-path")
solver <- policy_option_value(args, "--solver", "gurobi")
dir.create(allocation_path, recursive = TRUE, showWarnings = FALSE)
solve_checked <- function(edges, demand, population, target_rate, label) {
  if (any(!is.finite(demand)) || !is.finite(target_rate)) stop("Undefined median demand: ", label)
  maximum <- sum(vapply(split(demand, edges$village_i), max, numeric(1)) * population)
  if (maximum + 1e-5 < target_rate * sum(population)) {
    saveRDS(list(status = "target_infeasible", maximum = maximum, target = target_rate * sum(population)),
            file.path(allocation_path, paste0(label, ".rds")))
    return(NA_real_)
  }
  fit <- policy_cost_solve(edges = edges, demand = demand, population = population,
    target_rate = target_rate, site_cost = 1, solver = solver, solver_threads = 1L,
    solver_seed = 0L, work_path = tempdir(),
    log_file = file.path(allocation_path, paste0(label, ".log")))
  if (!isTRUE(fit$diagnostics$optimal)) stop("Uncertified median solution: ", label)
  saveRDS(fit, file.path(allocation_path, paste0(label, ".rds")))
  if (!is.null(canonical_policy_path) && label %in% c("cap-3500-control", "cap-3500-bracelet")) {
    directory <- file.path(canonical_policy_path, "median-allocations", sub("cap-3500-", "", label))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    saveRDS(fit, file.path(directory, "replicate-0001.rds"))
  }
  fit$summary$sites
}

if (!file.exists(compact_path)) {
  stop("Missing compact baseline draws: ", compact_path, call. = FALSE)
}

draws <- read.csv(compact_path, check.names = FALSE, stringsAsFactors = FALSE)
parameter_columns <- names(draws)[vapply(draws, is.numeric, logical(1))]
parameter_columns <- setdiff(
  parameter_columns,
  c("chain", "iteration", "draw", "replicate")
)
median_draw <- as.data.frame(as.list(vapply(
  draws[parameter_columns], median, numeric(1), na.rm = TRUE
)), check.names = FALSE)
parameter <- if (parameter_type == "canonical") {
  as.list(median_draw)
} else if (parameter_type == "raw") {
  as.list(canonical_policy_parameters(
    median_draw, replicate = 1L, mode_csv = compact_path
  ))
} else stop("--parameter-type must be raw or canonical.", call. = FALSE)
parameter$model_family <- "gaussian"
parameter$draw <- 1L
parameter$replicate <- 1L

distance_object <- readRDS(distance_path)
if (!identical(distance_object$candidate_site_mode, "all") ||
    nrow(distance_object$pot_df) != 1451L ||
    length(unique(distance_object$pot_df$cluster.id)) != 1451L) {
  stop("Distance-cap table requires the canonical all-1,451-site object.")
}
parameter$sd_of_dist <- distance_object$sd_of_dist
villages <- distance_object$village_df

population_table <- policy_adult_population(distance_object)
population <- population_table$population


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

  geographic_floor <- solve_checked(
    edges = edges, demand = rep(1, nrow(edges)),
    population = rep(1, nrow(villages)), target_rate = 0,
    label = paste0("cap-", cap, "-geographic-floor")
  )

  sites <- vapply(c("control", "bracelet"), function(regime) {
    demand <- predict_policy_draw(
      parameter, edges$distance,
      policy_scenarios[policy_scenarios$scenario == regime, , drop = FALSE]
    )$demand
    solve_checked(
      edges = edges, demand = demand, population = population,
      target_rate = target_rate, label = paste0("cap-", cap, "-", regime)
    )
  }, numeric(1))

  data.frame(
    cap_m = cap,
    feasible_links = nrow(edges),
    mean_options_per_village = mean(village_degree),
    villages_with_at_most_two_options = sum(village_degree <= 2),
    sites_feasible_for_multiple_villages = sum(site_degree >= 2),
    geographic_minimum_sites = geographic_floor,
    population_total = sum(population),
    target_expected_adults = target_rate * sum(population),
    control_status = if (is.na(sites["control"])) "target_infeasible" else "complete",
    bracelet_status = if (is.na(sites["bracelet"])) "target_infeasible" else "complete",
    control_sites = sites["control"],
    bracelet_sites = sites["bracelet"],
    sites_saved = sites["control"] - sites["bracelet"],
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, rows)

dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
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
