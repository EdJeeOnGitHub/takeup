#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")

input_path <- policy_option_value(args, "--input-path")
if (is.null(input_path)) stop("--input-path is required.", call. = FALSE)
parameter_status <- read.csv(
  file.path(input_path, "policy-model-parameter-status.csv"),
  stringsAsFactors = FALSE
)

scenario_results <- do.call(rbind, lapply(policy_scenarios$scenario, function(scenario) {
  path <- file.path(input_path, "allocations", scenario, "status.csv")
  if (!file.exists(path)) stop("Missing scenario status: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE)
}))
expected_draws <- parameter_status$draws[1]
if (any(table(scenario_results$scenario) != expected_draws)) {
  stop("One or more policy scenarios are incomplete.", call. = FALSE)
}

experimental_demand <- readRDS(file.path(input_path, "policy-experimental-demand.rds"))
experimental <- aggregate(
  cbind(mean_demand = experimental_demand$demand, mean_distance = experimental_demand$distance),
  by = list(draw = experimental_demand$draw, replicate = experimental_demand$replicate),
  FUN = mean
)
experimental$scenario_id <- 0L
experimental$scenario <- "experimental"
experimental$scenario_label <- "Experimental allocation"
experimental$status <- "observed_allocation"
experimental$solver_status <- NA_integer_
experimental$elapsed_seconds <- NA_real_
experimental$n_pot <- 144L
experimental$achieved_welfare <- experimental$mean_demand * 144
experimental$target_welfare <- unique(scenario_results$target_welfare)[1]
for (column in setdiff(names(scenario_results), names(experimental))) experimental[[column]] <- NA
all_results <- rbind(experimental[, names(scenario_results), drop = FALSE], scenario_results)
all_results$model_id <- parameter_status$model_id[1]
all_results$model_label <- parameter_status$model_label[1]
write.csv(all_results, file.path(input_path, "policy-model-replicates.csv"), row.names = FALSE)

summary_rows <- do.call(rbind, lapply(split(all_results, all_results$scenario_id), function(value) {
  n_pot <- summarize_policy_values(value$n_pot)
  takeup <- summarize_policy_values(value$mean_demand)
  distance <- summarize_policy_values(value$mean_distance / 1000)
  data.frame(
    model_id = value$model_id[1], model_label = value$model_label[1],
    scenario_id = value$scenario_id[1], scenario = value$scenario[1],
    scenario_label = value$scenario_label[1],
    n_pot_estimate = n_pot["estimate"], n_pot_low = n_pot["conf_low"],
    n_pot_high = n_pot["conf_high"], takeup_estimate = takeup["estimate"],
    takeup_low = takeup["conf_low"], takeup_high = takeup["conf_high"],
    distance_estimate = distance["estimate"], distance_low = distance["conf_low"],
    distance_high = distance["conf_high"], draws = nrow(value),
    target_infeasible_share = mean(value$status == "target_infeasible"),
    equilibrium_undefined_share = mean(value$status == "equilibrium_undefined"),
    stringsAsFactors = FALSE
  )
}))
summary_rows <- summary_rows[order(summary_rows$scenario_id), ]
write.csv(summary_rows, file.path(input_path, "policy-model-summary.csv"), row.names = FALSE)

paired_contrast <- function(left_scenario, right_scenario, label) {
  left <- all_results[all_results$scenario == left_scenario, ]
  right <- all_results[all_results$scenario == right_scenario, ]
  right <- right[match(left$draw, right$draw), ]
  if (any(left$draw != right$draw)) stop("Policy contrast draws do not align.", call. = FALSE)
  saved <- left$n_pot - right$n_pot
  interval <- summarize_policy_values(saved)
  data.frame(
    model_id = parameter_status$model_id[1],
    model_label = parameter_status$model_label[1],
    contrast = label,
    left_scenario = left_scenario,
    right_scenario = right_scenario,
    left_n_pot = median(left$n_pot, na.rm = TRUE),
    right_n_pot = median(right$n_pot, na.rm = TRUE),
    n_pot_saved = interval["estimate"],
    n_pot_saved_low = interval["conf_low"],
    n_pot_saved_high = interval["conf_high"],
    defined_draws = sum(is.finite(saved)),
    draws = length(saved),
    equilibrium_undefined_share = mean(
      left$status == "equilibrium_undefined" |
        right$status == "equilibrium_undefined"
    ),
    target_infeasible_share = mean(
      left$status == "target_infeasible" |
        right$status == "target_infeasible"
    ),
    stringsAsFactors = FALSE
  )
}
contrast_rows <- rbind(
  paired_contrast("control", "bracelet", "Endogenous social image returns"),
  paired_contrast("static-control", "static-bracelet", "Social image returns fixed at 0.5km")
)
write.csv(
  contrast_rows, file.path(input_path, "policy-model-contrast-summary.csv"),
  row.names = FALSE
)
message("Summarized ", parameter_status$model_label[1], " across ", expected_draws, " draws.")
