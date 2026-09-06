#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
source("R/policy/population.R")
args <- commandArgs(TRUE)
input <- policy_option_value(args, "--input-path")
root <- policy_option_value(args, "--shard-root")
map <- read.csv(file.path(input, "policy-edge-demand-draw-map.csv"))
edges <- readRDS(file.path(input, "policy-feasible-edges.rds"))
key <- paste(edges$village_i, edges$pot_j, sep = ":")
markers <- list.files(root, pattern = "^_SUCCESS$", recursive = TRUE, full.names = TRUE)
if (!length(markers)) stop("No successful shards.")
statuses <- list()
for (marker in markers) {
  shard <- dirname(marker)
  files <- list.files(file.path(shard, "allocations"), pattern = "^status.csv$", recursive = TRUE, full.names = TRUE)
  for (path in files) {
    status <- read.csv(path)
    scenario <- unique(status$scenario)
    if (length(scenario) != 1L || !scenario %in% policy_scenarios$scenario) stop("Invalid shard scenario.")
    dest <- file.path(input, "allocations", scenario)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_len(nrow(status))) {
      from <- file.path(dirname(path), sprintf("replicate-%04d.rds", status$draw[i]))
      result <- readRDS(from); a <- result$allocation; s <- result$status
      if (s$draw != status$draw[i] || s$status != status$status[i]) stop("Shard status disagrees with saved draw.")
      if (s$status != "equilibrium_undefined") {
        index <- match(paste(a$village_i, a$pot_j, sep = ":"), key)
        if (anyNA(index) || nrow(a) != 144L || anyDuplicated(a$village_i) ||
            any(a$distance != edges$distance[index]) ||
            abs(sum(a$expected_takers) - s$achieved_welfare) > 1e-7) stop("Invalid saved assignment.")
        if (s$status == "complete" && (!isTRUE(result$solver$optimal) ||
            sum(a$expected_takers) + 1e-5 < s$target_welfare)) stop("Uncertified or target-missing solution.")
      }
      target <- file.path(dest, basename(from))
      if (file.exists(target)) {
        if (tools::md5sum(target) != tools::md5sum(from)) stop("Conflicting assignment already collected.")
      } else if (!file.link(from, target) && !file.copy(from, target)) stop("Failed to publish assignment.")
    }
    statuses[[length(statuses) + 1L]] <- status
  }
}
all <- do.call(rbind, statuses)
if (anyDuplicated(all[c("scenario", "draw")]) || nrow(all) != 5L * nrow(map) ||
    !setequal(all$scenario, policy_scenarios$scenario) ||
    any(!all$status %in% c("complete", "target_infeasible", "equilibrium_undefined"))) stop("Incomplete/duplicate collection.")
for (scenario in policy_scenarios$scenario) {
  status <- all[all$scenario == scenario, ]; status <- status[order(status$draw), ]
  if (!identical(status$draw, sort(map$draw))) stop("Collected draw IDs differ from prediction inventory.")
  write.csv(status, file.path(input, "allocations", scenario, "status.csv"), row.names = FALSE)
}
write.csv(data.frame(draws = nrow(map), scenarios = 5L, assignments_checked = nrow(all), status = "complete"),
          file.path(input, "policy-collection-audit.csv"), row.names = FALSE)
