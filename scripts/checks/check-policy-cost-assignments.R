#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
root <- policy_option_value(args, "--run-root")
analysis <- policy_option_value(args, "--analysis-id")
report <- file.path(root, "reports")
model <- if (analysis == "baseline-posterior") "benchmark" else "cluster-weighted"
input <- file.path(root, "cap-3500", model)
distance <- readRDS(policy_option_value(args, "--distance-data"))
env <- new.env(); load(policy_option_value(args, "--census-data"), env)
census <- env$census.data
adults <- census[which(census$age.census >= 18 & census$cluster.id %in% distance$village_df$cluster.id), ]
stopifnot(nrow(adults) == 39301L, !anyDuplicated(adults$KEY.individ))
population <- as.numeric(table(factor(adults$cluster.id, levels = distance$village_df$cluster.id)))
edges <- distance$long_distance_mat
key <- paste(edges$index_i, edges$index_j)
targets <- read.csv(file.path(input, "policy-experimental-targets.csv"))
x <- read.csv(file.path(report, paste0("policy-allocation-draws-", analysis, ".csv")))
expected <- if (model == "benchmark") 1600L else 999L
stopifnot(length(unique(x$draw)) == expected, !anyDuplicated(x[c("draw", "estimand", "regime")]))
near <- function(a,b) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=1e-7,scale=1)))
paths <- character(nrow(x)); checked <- new.env(parent=emptyenv())
for (i in seq_len(nrow(x))) {
  row <- x[i, ]; legacy <- row$estimand == "legacy"
  weights <- if (legacy) rep(1,144) else population
  target <- targets[[if (legacy) "target_community_welfare" else "target_expected_adults"]][match(row$draw,targets$draw)]
  cache <- if (legacy) paste0("legacy-",row$regime) else if (row$regime == "bracelet" && is.finite(row$pooling_rho) && row$pooling_rho > 0) sprintf("population-bracelet-%03d",round(100*row$pooling_rho)) else NULL
  path <- if (is.null(cache)) file.path(input,"allocations",row$regime,sprintf("replicate-%04d.rds",row$draw)) else file.path(report,"sensitivity-allocations",analysis,cache,sprintf("replicate-%04d.rds",row$draw))
  saved <- if (exists(path,checked,inherits=FALSE)) get(path,checked) else readRDS(path)
  # Cache only within a draw to avoid retaining all allocation data in memory.
  if (i > 1 && x$draw[i] != x$draw[i-1]) rm(list=ls(checked),envir=checked)
  assign(path,saved,checked);paths[i] <- path
  a <- saved$allocation; solver <- if (is.null(cache)) saved$solver else saved$diagnostics
  stopifnot(nrow(a)==144L,!anyDuplicated(a$village_i),setequal(a$village_i,1:144),all(is.finite(a$demand)),isTRUE(solver$optimal))
  index <- match(paste(a$village_i,a$pot_j),key)
  stopifnot(!anyNA(index),all(edges$dist[index]<=3500))
  near(a$distance,edges$dist[index]);near(a$distance_km,edges$dist_km[index]);near(a$population,weights[a$village_i])
  takers <- weights[a$village_i]*a$demand
  near(a$expected_takers,takers);near(row$expected_takers,sum(takers))
  near(row$target_slack_takers,sum(takers)-target)
  stopifnot(sum(takers)+1e-4>=target,row$sites==length(unique(a$pot_j)),row$sites-solver$bound<1-1e-5)
  near(row$roundtrip_participant_km,sum(takers*2*a$distance_km))
  near(row$population_mean_distance_km,weighted.mean(a$distance_km,weights[a$village_i]))
}
costs <- read.csv(file.path(report,paste0("policy-break-even-draws-",analysis,".csv")))
control <- x[x$estimand=="population" & x$regime=="control", ];bracelet <- x[x$estimand=="population" & x$regime=="bracelet", ]
a <- control[match(costs$draw,control$draw), ];b <- bracelet[match(costs$draw,bracelet$draw), ]
near(costs$signal_cost_difference,.20*b$expected_takers)
near(costs$roundtrip_participant_km_difference,b$roundtrip_participant_km-a$roundtrip_participant_km)
for (price in c(100,250,500,1000)) near(costs[[paste0("cost_saved_at_",price)]],price*(a$sites-b$sites)-.20*b$expected_takers-costs$travel_cost*(b$roundtrip_participant_km-a$roundtrip_participant_km))
write.csv(data.frame(analysis_id=analysis,draws=expected,rows_checked=nrow(x),unique_assignments=length(unique(paths)),status="passed"),file.path(report,paste0("independent-cost-assignment-audit-",analysis,".csv")),row.names=FALSE)
message("PASS: ",analysis," assignments, targets, population weights, and cost grids.")
