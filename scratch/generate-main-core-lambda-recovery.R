#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(purrr)
  library(rlang)
})

workspace <- option_value("--workspace")
model_name <- option_value(
  "--model", "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)
fit_csvs <- strsplit(option_value("--fit-csvs", ""), ",", fixed = TRUE)[[1L]]
output_path <- option_value(
  "--output-path", "temp-data/main-core-lambda-identification/recovery"
)
stan_path <- option_value("--stan-path", "stan_models")
replicates <- as.integer(option_value("--replicates", "50"))
seed <- as.integer(option_value("--seed", "20263000"))
cmdstan_path_option <- option_value(
  "--cmdstan-path", Sys.getenv("CMDSTAN", unset = "")
)
if (is.null(workspace) || !file.exists(workspace) ||
    length(fit_csvs) < 2L || any(!file.exists(fit_csvs))) {
  stop("Existing workspace and baseline fit CSVs are required.", call. = FALSE)
}
if (nzchar(cmdstan_path_option)) {
  Sys.setenv(CMDSTAN = cmdstan_path_option)
  set_cmdstan_path(cmdstan_path_option)
}
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE")
options(cmdstanr_no_ver_check = TRUE)

data <- prepare_main_core_data(
  workspace,
  model_name,
  lambda_structure = 0L,
  lambda_log_ratio_sd_prior = 0.25
)
fit <- read_cmdstan_csv(fit_csvs)
draws <- as_draws_df(fit$post_warmup_draws)

# Select the coherent posterior draw nearest the componentwise median of the
# parameters that matter for the structural decomposition. Unlike a vector of
# componentwise medians, an observed joint draw preserves simplex constraints.
medoid_variables <- intersect(
  c(
    "base_mu_rep", "beta_intercept", "beta_ink_effect",
    "beta_bracelet_effect", "wtp_value_utility", "hyper_wtp_mu",
    "wtp_sigma", "dist_beta_v[1]", "raw_u_sd[1]",
    paste0("hyper_beta_1ord[", 1:4, "]"),
    paste0("hyper_dist_beta_1ord[", 1:4, "]")
  ),
  names(draws)
)
medoid_matrix <- as.matrix(draws[, medoid_variables, drop = FALSE])
medoid_center <- apply(medoid_matrix, 2, median)
medoid_scale <- apply(medoid_matrix, 2, mad)
medoid_scale[!is.finite(medoid_scale) | medoid_scale == 0] <- 1
medoid_distance <- rowSums(sweep(
  sweep(medoid_matrix, 2, medoid_center, "-"),
  2,
  medoid_scale,
  "/"
)^2)
medoid_index <- which.min(medoid_distance)

# Construct a valid one-draw CmdStan CSV by copying the exact raw row at that
# medoid. Raw CmdStan headers use dotted array indices whereas posterior uses
# bracketed indices, so reconstructing the row by column name is unsafe.
medoid_chain <- draws$.chain[medoid_index]
medoid_iteration <- draws$.iteration[medoid_index]
template_lines <- readLines(fit_csvs[medoid_chain])
header_index <- which(!startsWith(template_lines, "#"))[1]
data_lines <- which(seq_along(template_lines) > header_index &
                    !startsWith(template_lines, "#") &
                    nzchar(template_lines))
if (medoid_iteration < 1L || medoid_iteration > length(data_lines)) {
  stop("Could not locate the posterior-medoid CSV row.", call. = FALSE)
}
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
median_csv <- file.path(output_path, "posterior-medoid-parameters.csv")
writeLines(
  c(
    template_lines[seq_len(header_index)],
    template_lines[data_lines[medoid_iteration]]
  ),
  median_csv
)

scenarios <- data.frame(
  scenario = c("common", "grouped-025", "grouped-050", "arm-heterogeneous", "no-image"),
  zero_lambda = c(0L, 0L, 0L, 0L, 1L),
  stringsAsFactors = FALSE
)
deviations <- list(
  c(0, 0, 0, 0),
  c(-0.125, 0.125, -0.125, 0.125),
  c(-0.25, 0.25, -0.25, 0.25),
  c(-0.10, -0.20, 0.10, 0.20),
  c(0, 0, 0, 0)
)
sim_model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core_simulate.stan"),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
compact_model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core_compact_gq.stan"),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
manifest <- list()
for (scenario_index in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[scenario_index]
  scenario_path <- file.path(output_path, "datasets", scenario)
  dir.create(scenario_path, recursive = TRUE, showWarnings = FALSE)
  sim_data <- data
  sim_data$core_sim_zero_lambda <- scenarios$zero_lambda[scenario_index]
  sim_data$core_sim_lambda_log_deviation <- deviations[[scenario_index]]
  true_lambda <- if (scenarios$zero_lambda[scenario_index]) {
    rep(0, 4)
  } else {
    draws$base_mu_rep[medoid_index] * exp(deviations[[scenario_index]])
  }
  true_data <- data
  true_data$core_gq_override_lambda <- 1L
  true_data$core_gq_lambda_override <- true_lambda
  true_gq <- compact_model$generate_quantities(
    fitted_params = median_csv,
    data = true_data,
    seed = seed + scenario_index,
    output_dir = scenario_path,
    output_basename = "true-estimands",
    threads_per_chain = 1
  )
  true_draw <- as_draws_df(true_gq$draws(variables = c(
    "core_compact_sm_rescaled", "core_compact_takeup_level"
  )))
  true_sm_1500 <- -as.numeric(unlist(true_draw[1, paste0(
    "core_compact_sm_rescaled[2,", 1:4, "]"
  )], use.names = FALSE))
  for (replicate in seq_len(replicates)) {
    replicate_seed <- seed + scenario_index * 10000L + replicate
    gq <- sim_model$generate_quantities(
      fitted_params = median_csv,
      data = sim_data,
      seed = replicate_seed,
      output_dir = scenario_path,
      output_basename = sprintf("sim-%03d", replicate),
      threads_per_chain = 1
    )
    sim_draw <- as_draws_df(gq$draws(variables = c(
      "core_sim_takeup", "core_sim_num_knows_1ord",
      "core_sim_num_knows_2ord", "core_sim_wtp_response",
      "core_sim_gift_choice", "core_sim_signal_lambda"
    )))
    extract_array <- function(prefix, length) {
      as.integer(unlist(
        sim_draw[1, paste0(prefix, "[", seq_len(length), "]")],
        use.names = FALSE
      ))
    }
    refit_data <- data
    refit_data$takeup <- extract_array("core_sim_takeup", data$num_obs)
    refit_data$num_knows_1ord <- extract_array(
      "core_sim_num_knows_1ord", data$num_beliefs_obs
    )
    refit_data$num_knows_2ord <- extract_array(
      "core_sim_num_knows_2ord", data$num_beliefs_obs
    )
    refit_data$wtp_response <- extract_array(
      "core_sim_wtp_response", data$num_wtp_obs
    )
    refit_data$gift_choice <- extract_array(
      "core_sim_gift_choice", data$num_wtp_obs
    )
    # Preserve Stan's vector dimension when this fit-105 setting contains only
    # one optimization distance. A bare length-one R vector is serialized as a
    # JSON scalar by write_stan_json().
    if (length(refit_data$optim_distances) == 1L &&
        is.null(dim(refit_data$optim_distances))) {
      refit_data$optim_distances <- array(refit_data$optim_distances, dim = 1L)
    }
    data_json <- file.path(scenario_path, sprintf("rep-%03d.json", replicate))
    write_stan_json(refit_data, data_json)
    manifest[[length(manifest) + 1L]] <- data.frame(
      scenario = scenario,
      replicate = replicate,
      seed = replicate_seed,
      data_json = data_json,
      true_lambda_control = true_lambda[1],
      true_lambda_ink = true_lambda[2],
      true_lambda_calendar = true_lambda[3],
      true_lambda_bracelet = true_lambda[4],
      true_sm_control_1500m = true_sm_1500[1],
      true_sm_ink_1500m = true_sm_1500[2],
      true_sm_calendar_1500m = true_sm_1500[3],
      true_sm_bracelet_1500m = true_sm_1500[4],
      true_signal_log_ratio = if (all(true_lambda > 0)) {
        mean(log(true_lambda[c(2, 4)])) -
          mean(log(true_lambda[c(1, 3)]))
      } else {
        NA_real_
      }
    )
  }
}
manifest <- bind_rows(manifest)
fit_structures <- data.frame(
  fit_structure = c("common", "grouped", "arm"),
  lambda_structure = 0:2,
  prior_sd = c(0.25, 0.25, 0.25)
)
fit_manifest <- merge(manifest, fit_structures, by = NULL)
fit_manifest$task_id <- seq_len(nrow(fit_manifest))
fit_manifest$label <- sprintf(
  "%s-rep%03d-%s",
  fit_manifest$scenario,
  fit_manifest$replicate,
  fit_manifest$fit_structure
)
fit_manifest$status <- "pending"
write.csv(
  fit_manifest,
  file.path(output_path, "lambda-recovery-manifest.csv"),
  row.names = FALSE
)
print(table(fit_manifest$scenario, fit_manifest$fit_structure))
