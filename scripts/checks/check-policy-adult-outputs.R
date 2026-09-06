#!/usr/bin/env Rscript
# Independent audit: recompute counts/targets/objectives from source observations.
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
root <- policy_option_value(args, "--run-root")
distance <- readRDS(policy_option_value(args, "--distance-data"))
stopifnot(identical(distance$candidate_site_mode, "all"), nrow(distance$pot_df) == 1451L,
          !anyDuplicated(distance$pot_df$cluster.id), nrow(distance$long_distance_mat) == 144L * 1451L,
          !anyDuplicated(distance$long_distance_mat[c("index_i", "index_j")]))
env <- new.env(); load(policy_option_value(args, "--census-data"), env)
census <- env$census.data
adults <- census[which(census$age.census >= 18 & census$cluster.id %in% distance$village_df$cluster.id), ]
stopifnot(nrow(adults) == 39301L, !anyDuplicated(adults$KEY.individ))
population <- as.numeric(table(factor(adults$cluster.id, levels = distance$village_df$cluster.id)))
stopifnot(length(population) == 144L, all(population > 0))
paths <- list.files(root, pattern = "^policy-population.csv$", recursive = TRUE, full.names = TRUE)
model_filter <- policy_option_value(args, "--model")
cap_filter <- policy_option_value(args, "--cap")
if (!is.null(model_filter)) paths <- paths[basename(dirname(paths)) == model_filter]
if (!is.null(cap_filter)) paths <- paths[basename(dirname(dirname(paths))) == paste0("cap-", cap_filter)]
output_path <- policy_option_value(args, "--output-path", root)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
expected <- as.integer(policy_option_value(args, "--expected-combinations", "13"))
draw_limit <- as.integer(policy_option_value(args, "--max-draws", "0"))
draw_counts <- c(benchmark = 1600L, `private-distance-community-image` = 800L,
  `full-information` = 800L, `exclude-dispersed` = 4000L, `cluster-shock` = 4000L,
  `tight-multinomial` = 1600L, `second-order-observability` = 4000L,
  `grouped-lambda` = 1600L, `arm-lambda` = 1600L, `student-t5` = 1600L,
  `finite-mixture` = 3200L, `cluster-weighted` = 999L)
stopifnot(length(paths) == expected)
near <- function(a, b, tolerance = 1e-7) {
  if (!isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance, scale = 1))) stop("Numeric audit mismatch")
}
audit <- list(); benchmark_targets <- list()
for (path in paths) {
  input <- dirname(path); message("Auditing ", input)
  cap <- as.numeric(sub("cap-([0-9]+).*", "\\1", basename(dirname(input))))
  p <- read.csv(path); stopifnot(identical(p$village_i, 1:144)); near(p$population, population)
  target <- read.csv(file.path(input, "policy-experimental-targets.csv"))
  experimental <- readRDS(file.path(input, "policy-experimental-demand.rds"))
  map <- read.csv(file.path(input, "policy-edge-demand-draw-map.csv"))
  parameters <- read.csv(file.path(input, "policy-model-parameters.csv"))
  comparison_path <- policy_option_value(args, "--comparison-200-csv")
  if (basename(input) == "benchmark" && !is.null(comparison_path)) {
    reference <- read.csv(comparison_path)
    index <- match(paste(reference$chain, reference$iteration), paste(parameters$chain, parameters$iteration))
    stopifnot(nrow(reference) == 200L, !anyNA(index), sum(parameters$comparison_200) == 200L,
              all(parameters$comparison_200[index]))
    columns <- setdiff(names(reference)[vapply(reference, is.numeric, logical(1))], c("draw", "replicate"))
    stopifnot(all(columns %in% names(parameters)))
    for (column in columns) near(reference[[column]], parameters[[column]][index], 1e-12)
  }
  count <- unname(draw_counts[basename(input)])
  stopifnot(length(count) == 1L, !is.na(count), nrow(parameters) == count,
            nrow(map) == if (draw_limit > 0L) min(draw_limit, count) else count,
            !anyDuplicated(map$draw), !anyDuplicated(parameters$draw))
  edges <- readRDS(file.path(input, "policy-feasible-edges.rds"))
  original <- distance$long_distance_mat[distance$long_distance_mat$dist <= cap, ]
  stopifnot(nrow(edges) == nrow(original), setequal(paste(edges$village_i, edges$pot_j), paste(original$index_i, original$index_j)))
  source_index <- match(paste(edges$village_i, edges$pot_j), paste(original$index_i, original$index_j))
  near(edges$distance, original$dist[source_index])
  demand <- readRDS(file.path(input, "policy-edge-demand-matrix.rds"))
  stopifnot(!anyDuplicated(target$draw), setequal(target$draw, map$draw), all(edges$distance <= cap))
  edge_keys <- paste(edges$village_i, edges$pot_j)
  for (draw in map$draw) {
    e <- experimental[experimental$draw == draw, ]; stopifnot(nrow(e) == 144L, !anyDuplicated(e$village_i))
    actual_target <- sum(e$demand * population[e$village_i])
    near(actual_target, target$target_expected_adults[target$draw == draw])
    for (sid in 1:5) {
      scenario <- policy_scenarios$scenario[sid]
      saved <- readRDS(file.path(input, "allocations", scenario, sprintf("replicate-%04d.rds", draw)))
      s <- saved$status; a <- saved$allocation
      stopifnot(s$draw == draw, s$scenario == scenario)
      near(s$target_welfare, actual_target)
      d <- demand[match(draw, map$draw), (sid - 1L) * nrow(edges) + seq_len(nrow(edges))]
      undefined <- any(!is.finite(d)) || !is.finite(actual_target)
      if (undefined) {
        stopifnot(s$status == "equilibrium_undefined", nrow(a) == 0L)
      } else {
        stopifnot(nrow(a) == 144L, setequal(a$village_i, 1:144), !anyDuplicated(a$village_i))
        index <- match(paste(a$village_i, a$pot_j), edge_keys); stopifnot(!anyNA(index))
        near(a$demand, d[index]); near(a$distance, edges$distance[index]); near(a$population, population[a$village_i])
        takers <- population[a$village_i] * d[index]
        near(a$expected_takers, takers); near(s$achieved_welfare, sum(takers))
        near(s$mean_distance, weighted.mean(a$distance, a$population))
        near(s$mean_demand, sum(takers) / 39301)
        stopifnot(s$n_pot == length(unique(a$pot_j)))
        maximum <- sum(vapply(split(d, edges$village_i), max, numeric(1)) * population)
        if (maximum + 1e-5 < actual_target) {
          stopifnot(s$status == "target_infeasible"); near(sum(takers), maximum)
        } else {
          stopifnot(s$status == "complete", isTRUE(saved$solver$optimal), sum(takers) + 1e-5 >= actual_target)
          # Gurobi's binary feasibility tolerance can leave sub-micro-unit
          # objective residue; the bound below separately proves the integer optimum.
          near(saved$solver$objective, s$n_pot, tolerance = 1e-5)
          stopifnot(is.finite(saved$solver$gap), saved$solver$gap <= 1e-4,
                    s$n_pot - saved$solver$bound < 1 - 1e-5)
        }
      }
      audit[[length(audit) + 1L]] <- data.frame(model = basename(input), cap = cap, draw = draw,
        scenario = scenario, status = s$status, sites = s$n_pot, seconds = s$elapsed_seconds)
    }
  }
  if (basename(input) == "benchmark") benchmark_targets[[as.character(cap)]] <- target[order(target$draw), c("draw", "target_expected_adults")]
}
if (length(benchmark_targets) > 1L) for (x in benchmark_targets[-1L]) stopifnot(identical(x, benchmark_targets[[1L]]))
write.csv(do.call(rbind, audit), file.path(output_path, "independent-allocation-audit.csv"), row.names = FALSE)
writeLines("PASS: census counts, draw targets, caps, demand joins, weighted summaries, infeasibility, solver certificates, cross-cap targets",
           file.path(output_path, "independent-allocation-audit.txt"))
