#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
root <- policy_option_value(args, "--run-root")
output <- policy_option_value(args, "--output-path", file.path(root, "reports"))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
models <- c("benchmark", "private-distance-community-image", "full-information", "exclude-dispersed",
            "cluster-shock", "tight-multinomial", "second-order-observability", "grouped-lambda",
            "arm-lambda", "student-t5", "finite-mixture", "cluster-weighted")
inventory <- rbind(data.frame(model = models, cap = 3500L), data.frame(model = "benchmark", cap = c(4500L,5500L,10000L)))
draws <- list(); paired <- list(); summaries <- list()
interval <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(low = NA_real_, median = NA_real_, high = NA_real_))
  setNames(as.numeric(quantile(x, c(.025,.5,.975))), c("low", "median", "high"))
}
for (i in seq_len(nrow(inventory))) {
  item <- inventory[i, ]; input <- file.path(root, paste0("cap-", item$cap), item$model)
  x <- read.csv(file.path(input, "policy-model-replicates.csv"))
  x$cap_m <- item$cap
  parameters <- read.csv(file.path(input, "policy-model-parameters.csv"))
  x$comparison_200 <- if ("comparison_200" %in% names(parameters)) parameters$comparison_200[match(x$draw, parameters$draw)] else FALSE
  stopifnot(!anyDuplicated(x[c("scenario", "draw")]), all(x$population_total == 39301),
            all(x$population_weighting == "adult-census"))
  draws[[i]] <- x
  control <- x[x$scenario == "control", ]
  for (scenario in policy_scenarios$scenario) {
    y <- x[x$scenario == scenario, ]; y <- y[match(control$draw, y$draw), ]
    stopifnot(identical(y$draw, control$draw))
    valid <- control$status == "complete" & y$status == "complete"
    delta <- control$n_pot - y$n_pot
    z <- data.frame(model_id = item$model, cap_m = item$cap, scenario = scenario,
      draw = y$draw, replicate = y$replicate, comparison_200 = y$comparison_200,
      control_status = control$status, scenario_status = y$status,
      both_targets_met = valid, sites_saved_including_infeasible_fallbacks = delta,
      sites_saved_at_preserved_target = ifelse(valid, delta, NA_real_),
      takeup_change = y$mean_demand - control$mean_demand,
      adult_mean_distance_change_km = (y$mean_distance - control$mean_distance) / 1000)
    paired[[length(paired) + 1L]] <- z
    q <- interval(z$sites_saved_at_preserved_target)
    sites <- interval(y$n_pot[y$status == "complete"])
    summaries[[length(summaries) + 1L]] <- data.frame(model_id = item$model, cap_m = item$cap,
      scenario = scenario, draws = nrow(y), target_met = sum(y$status == "complete"),
      target_infeasible = sum(y$status == "target_infeasible"),
      equilibrium_undefined = sum(y$status == "equilibrium_undefined"), paired_targets_met = sum(valid),
      sites_low = sites[1], sites_median = sites[2], sites_high = sites[3],
      paired_sites_saved_low = q[1], paired_sites_saved_median = q[2], paired_sites_saved_high = q[3],
      probability_sites_saved_given_both_targets_met = if (any(valid)) mean(delta[valid] > 0) else NA_real_)
  }
}
all <- do.call(rbind, draws)
stopifnot(nrow(do.call(rbind, summaries)) == 75L)
write.csv(all, file.path(output, "all-model-cap-replicates.csv"), row.names = FALSE)
write.csv(do.call(rbind, paired), file.path(output, "all-model-cap-paired-contrasts.csv"), row.names = FALSE)
write.csv(do.call(rbind, summaries), file.path(output, "all-model-cap-summary.csv"), row.names = FALSE)
benchmark <- all[all$model_id == "benchmark", ]
baseline <- benchmark[benchmark$cap_m == 3500, ]
caps <- do.call(rbind, lapply(c(4500,5500,10000), function(cap) {
  wider <- benchmark[benchmark$cap_m == cap, ]
  x <- merge(baseline, wider, by = c("draw", "replicate", "scenario"), suffixes = c("_3500", "_wider"))
  stopifnot(nrow(x) == nrow(baseline), all(abs(x$target_welfare_3500 - x$target_welfare_wider) < 1e-7))
  data.frame(draw = x$draw, scenario = x$scenario, cap_m = cap,
    status_3500 = x$status_3500, status_wider = x$status_wider,
    sites_saved = x$n_pot_3500 - x$n_pot_wider,
    distance_change_km = (x$mean_distance_wider - x$mean_distance_3500)/1000)
}))
write.csv(caps, file.path(output, "benchmark-paired-distance-cap-contrasts.csv"), row.names = FALSE)
write.csv(benchmark[benchmark$comparison_200, ], file.path(output, "benchmark-original-200-subset.csv"), row.names = FALSE)
