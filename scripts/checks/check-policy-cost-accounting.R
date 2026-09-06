#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
reference <- policy_option_value(args, "--reference-commit", "8486c9c3bd416f6b04902c705f0c075851d33564")
output <- policy_option_value(args, "--output", "temp-data/policy-adult-rerun/accounting-equivalence.csv")
path <- tempfile(fileext = ".R")
stopifnot(system2("git", c("show", paste0(reference, ":R/policy/cost-sensitivity.R")), stdout = path) == 0L)
old <- new.env(); sys.source(path, old); unlink(path)
new <- new.env(); sys.source("R/policy/cost-sensitivity.R", new)
edges <- expand.grid(village_i = 1:3, pot_j = 1:3)
edges$distance_km <- c(0,1,2,1,0,1,2,1,0)
edges$distance <- edges$distance_km * 1000
results <- lapply(c(FALSE, TRUE), function(pooling) {
  demand <- 0.9 - edges$distance_km * 0.1
  options <- list(edges = edges, demand = demand, population = c(10,20,30), target_rate = .7,
    site_cost = 2, signal_cost_per_taker = .03, travel_cost_per_roundtrip_km = .01,
    solver = "glpk", work_path = tempdir(), pooled_demand = if (pooling) demand * .95 else NULL)
  before <- do.call(old$policy_cost_solve, options)
  after <- do.call(new$policy_cost_solve, options)
  data.frame(reference = reference, pooling = pooling,
    allocations_identical = identical(before$allocation, after$allocation),
    summaries_identical = identical(before$summary, after$summary))
})
results <- do.call(rbind, results)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output, row.names = FALSE)
stopifnot(all(results$allocations_identical), all(results$summaries_identical))
message("PASS: both accounting regression cases preserve allocations and all summaries.")
