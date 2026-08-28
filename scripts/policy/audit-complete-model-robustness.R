#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")
package_path <- policy_option_value(args, "--package-path")
if (is.null(package_path) || !dir.exists(package_path)) {
  stop("--package-path must name an assembled package.", call. = FALSE)
}

catalog <- data.frame(
  model_id = c(
    "benchmark", "private-distance-community-image", "full-information",
    "exclude-dispersed", "cluster-shock", "tight-multinomial",
    "second-order-observability", "grouped-lambda", "arm-lambda",
    "student-t5", "finite-mixture"
  ),
  reader_facing_label = c(
    "Benchmark", "Individual travel costs", "Individual distance observed by peers",
    "Excluding geographically dispersed communities", "Unobserved community heterogeneity",
    "Correct classification of take-up", "Perceived community observability",
    "By public-signal status", "By treatment arm", "Heavy-tailed v distribution",
    "Mixture v distribution"
  ),
  source_fit_directory = c(
    "200 balanced draws from assigned-distance slim chains",
    "/project/akaring/takeup-data/data/stan_analysis_data/streamlined-active-robustness/private-distance-community-image/fits",
    "/project/akaring/takeup-data/data/stan_analysis_data/streamlined-active-robustness/full-information/fits",
    "/project/akaring/takeup-data/data/stan_analysis_data/streamlined-active-robustness/exclude-dispersed/fits",
    "/project/akaring/takeup-data/data/stan_analysis_data/main-core-cluster-shock-production",
    "/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors/tight",
    "/project/akaring/takeup-data/data/stan_analysis_data/streamlined-active-robustness/second-order-observability/fits",
    "/project/akaring/takeup-data/data/stan_analysis_data/main-core-lambda-identification/fits/grouped-sd0p25",
    "/project/akaring/takeup-data/data/stan_analysis_data/main-core-lambda-identification/fits/arm-sd0p25",
    "/project/akaring/takeup-data/data/stan_analysis_data/main-core-student-t-robustness/fits/student-t5",
    "/project/akaring/takeup-data/data/stan_analysis_data/main-core-finite-mixture-robustness-800/fits/finite-mixture"
  ), stringsAsFactors = FALSE
)

model_rows <- list()
scenario_rows <- list()
for (index in seq_len(nrow(catalog))) {
  item <- catalog[index, ]
  directory <- file.path(package_path, item$model_id)
  required <- c(
    "policy-model-summary.csv", "policy-model-contrast-summary.csv",
    "policy-model-parameter-status.csv", "policy-prediction-status.csv",
    "policy-model-parameters.csv", "allocations"
  )
  missing <- required[!file.exists(file.path(directory, required)) &
    !dir.exists(file.path(directory, required))]
  if (length(missing)) {
    stop(item$model_id, " missing: ", paste(missing, collapse = ", "))
  }
  parameter <- read.csv(file.path(directory, "policy-model-parameter-status.csv"))
  prediction <- read.csv(file.path(directory, "policy-prediction-status.csv"))
  summary <- read.csv(file.path(directory, "policy-model-summary.csv"))
  statuses <- do.call(rbind, lapply(policy_scenarios$scenario, function(scenario) {
    path <- file.path(directory, "allocations", scenario, "status.csv")
    if (!file.exists(path)) stop("Missing status file: ", path)
    value <- read.csv(path, stringsAsFactors = FALSE)
    value$model_id <- item$model_id
    value
  }))
  scenario_audit <- do.call(rbind, lapply(split(statuses, statuses$scenario), function(value) {
    data.frame(
      model_id = item$model_id, scenario = value$scenario[1],
      attempted_draws = nrow(value), expected_draws = parameter$draws[1],
      complete = nrow(value) == parameter$draws[1],
      undefined = sum(value$status == "equilibrium_undefined"),
      undefined_share = mean(value$status == "equilibrium_undefined"),
      target_infeasible = sum(value$status == "target_infeasible"),
      target_infeasible_share = mean(value$status == "target_infeasible"),
      stringsAsFactors = FALSE
    )
  }))
  scenario_rows[[index]] <- scenario_audit
  note <- switch(item$model_id,
    "private-distance-community-image" = paste(
      "Complete 9805-observation fit104 coordinate frame; community fixed point",
      "with household-specific travel-cost adjustment; 3.5km centroid assignment cap."
    ),
    "full-information" = paste(
      "Complete 9805-observation fit104 coordinate frame; household-specific",
      "distance and fixed point; five communities transported out of fit105 estimation sample."
    ),
    "exclude-dispersed" = "Alternative fit transported to the common 144-community policy geography.",
    "cluster-shock" = paste(
      "Reused validated 4000-draw cluster-shock parameters; prediction and",
      "optimization rerun on the common 1451-site geography."
    ),
    "benchmark" = "Exactly 200 balanced posterior draws; not the 210-refit bootstrap.",
    ""
  )
  model_rows[[index]] <- data.frame(
    model_id = item$model_id, reader_facing_label = item$reader_facing_label,
    source_fit_directory = item$source_fit_directory, refit_performed = "no",
    chains = parameter$chains[1], posterior_draws = parameter$draws[1],
    prediction_complete = nrow(summary) == 6L,
    optimization_complete = all(scenario_audit$complete),
    number_undefined = sum(scenario_audit$undefined),
    share_undefined = sum(scenario_audit$undefined) / nrow(statuses),
    number_target_infeasible = sum(scenario_audit$target_infeasible),
    share_target_infeasible = sum(scenario_audit$target_infeasible) / nrow(statuses),
    notes = note, stringsAsFactors = FALSE
  )
  provenance_path <- file.path(directory, "provenance.txt")
  if (!file.exists(provenance_path)) {
    writeLines(c(
      paste0("model_id=", item$model_id),
      paste0("reader_facing_label=", item$reader_facing_label),
      paste0("source_fit_directory=", item$source_fit_directory),
      "distance_definition=assigned", "structural_refit_performed=no",
      "candidate_sites=1451", "distance_cap_m=3500",
      paste0("posterior_draws=", parameter$draws[1]),
      paste0("chains=", parameter$chains[1])
    ), provenance_path)
  }
}

audit_path <- file.path(package_path, "audit")
dir.create(audit_path, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, model_rows), file.path(audit_path, "model-status.csv"), row.names = FALSE)
write.csv(do.call(rbind, scenario_rows), file.path(audit_path, "scenario-status.csv"), row.names = FALSE)
writeLines(c(
  "distance_definition=assigned", "structural_refits_performed=none",
  "candidate_sites=1451", "distance_cap_m=3500",
  paste0("git_commit=", system("git rev-parse HEAD", intern = TRUE)),
  paste0("generated_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(audit_path, "provenance.txt"))
message("Audited all eleven policy models in ", package_path)
