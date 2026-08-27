#!/usr/bin/env Rscript

# Prepare a conditional-on-design Monte Carlo for the assigned-distance
# benchmark model. This script simulates outcomes but never refits a model.

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(purrr)
  library(rlang)
})

workspace <- option_value(
  "--workspace", "data/stan_analysis_data/dist_fit104.RData"
)
fit_csvs <- strsplit(option_value("--fit-csvs", ""), ",", fixed = TRUE)[[1L]]
output_path <- option_value(
  "--output-path", "temp-data/main-core-benchmark-recovery"
)
stan_path <- option_value("--stan-path", "stan_models")
replicates <- as.integer(option_value("--replicates", "50"))
grid_points <- as.integer(option_value("--grid-points", "21"))
seed <- as.integer(option_value("--seed", "20260827"))
cmdstan_path_option <- option_value(
  "--cmdstan-path", Sys.getenv("CMDSTAN_PATH", unset = "")
)

if (!file.exists(workspace)) stop("Workspace not found: ", workspace)
if (length(fit_csvs) != 4L || any(!file.exists(fit_csvs))) {
  stop("--fit-csvs must list the four assigned-distance benchmark chains.")
}
if (replicates < 1L || grid_points < 5L || grid_points %% 2L != 1L) {
  stop("Use at least one replicate and an odd likelihood grid of at least five points.")
}
if (nzchar(cmdstan_path_option)) set_cmdstan_path(cmdstan_path_option)
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE")
options(cmdstanr_no_ver_check = TRUE)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dataset_path <- file.path(output_path, "datasets")
dir.create(dataset_path, recursive = TRUE, showWarnings = FALSE)

data <- prepare_main_core_data(
  workspace_path = workspace,
  model_name = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP",
  distance_definition = "assigned"
)
fit <- read_cmdstan_csv(fit_csvs)
draws <- as_draws_df(fit$post_warmup_draws)

# Choose one observed joint draw close to the posterior center. Keeping an
# actual row preserves all constraints and correlations among nuisance terms.
medoid_variables <- intersect(c(
  "base_mu_rep", "beta_intercept", "beta_ink_effect",
  "beta_bracelet_effect", "wtp_value_utility", "hyper_wtp_mu",
  "wtp_sigma", "dist_beta_v[1]", "raw_u_sd[1]",
  paste0("hyper_beta_1ord[", 1:4, "]"),
  paste0("hyper_dist_beta_1ord[", 1:4, "]")
), names(draws))
medoid_matrix <- as.matrix(draws[, medoid_variables, drop = FALSE])
center <- apply(medoid_matrix, 2, median)
scale <- apply(medoid_matrix, 2, mad)
scale[!is.finite(scale) | scale == 0] <- 1
distance <- rowSums(sweep(sweep(medoid_matrix, 2, center, "-"), 2, scale, "/")^2)
medoid_index <- which.min(distance)
medoid_chain <- draws$.chain[medoid_index]
medoid_iteration <- draws$.iteration[medoid_index]

template <- readLines(fit_csvs[medoid_chain])
header_index <- which(!startsWith(template, "#"))[1]
data_lines <- which(seq_along(template) > header_index &
                      !startsWith(template, "#") & nzchar(template))
if (medoid_iteration > length(data_lines)) stop("Cannot locate medoid CSV row.")
truth_row <- template[data_lines[medoid_iteration]]
truth_csv <- file.path(output_path, "truth-parameters.csv")
writeLines(c(template[seq_len(header_index)], truth_row), truth_csv)

# Build one multi-row fitted-parameter CSV containing identified normalized
# likelihood slices plus the raw-scale invariance path. CmdStan can evaluate
# this complete diagnostic grid in one GQ invocation per simulated dataset.
raw_header <- strsplit(template[header_index], ",", fixed = TRUE)[[1L]]
raw_truth <- strsplit(truth_row, ",", fixed = TRUE)[[1L]]
required_columns <- c(
  "dist_beta_v.1", "base_mu_rep", "raw_u_sd.1", "wtp_value_utility",
  "beta_intercept", "beta_ink_effect", "beta_calendar_effect",
  "beta_bracelet_effect"
)
if (!all(required_columns %in% raw_header)) stop("Benchmark CSV schema drift.")
truth_value <- function(parameter) {
  as.numeric(raw_truth[match(parameter, raw_header)])
}
set_value <- function(row, parameter, value) {
  row[match(parameter, raw_header)] <-
    format(value, scientific = TRUE, digits = 17)
  row
}
factors <- seq(0.5, 1.5, length.out = grid_points)
normalized_grid <- expand.grid(
  object = c("normalized_delta", "normalized_lambda"), factor = factors,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
s_truth <- sqrt(1 + truth_value("raw_u_sd.1")^2)
minimum_k <- max(0.75, 1 / s_truth + 0.01)
scale_factors <- sort(unique(c(seq(minimum_k, 1.5, length.out = grid_points), 1)))
scale_grid <- data.frame(
  object = "raw_scale_path", factor = scale_factors,
  stringsAsFactors = FALSE
)
grid <- bind_rows(normalized_grid, scale_grid)
grid$candidate_id <- seq_len(nrow(grid))
grid$truth <- ifelse(
  grid$object == "normalized_delta", truth_value("dist_beta_v.1") / s_truth,
  ifelse(grid$object == "normalized_lambda",
         truth_value("base_mu_rep") / s_truth^2, 1)
)
grid$value <- ifelse(
  grid$object == "raw_scale_path", grid$factor, grid$truth * grid$factor
)
candidate_rows <- vapply(seq_len(nrow(grid)), function(index) {
  row <- raw_truth
  factor <- grid$factor[index]
  if (grid$object[index] == "normalized_delta") {
    row <- set_value(row, "dist_beta_v.1", truth_value("dist_beta_v.1") * factor)
  } else if (grid$object[index] == "normalized_lambda") {
    row <- set_value(row, "base_mu_rep", truth_value("base_mu_rep") * factor)
  } else {
    row <- set_value(row, "raw_u_sd.1", sqrt((factor * s_truth)^2 - 1))
    row <- set_value(row, "dist_beta_v.1", truth_value("dist_beta_v.1") * factor)
    row <- set_value(row, "base_mu_rep", truth_value("base_mu_rep") * factor^2)
    for (parameter in c(
      "wtp_value_utility", "beta_intercept", "beta_ink_effect",
      "beta_calendar_effect", "beta_bracelet_effect"
    )) {
      row <- set_value(row, parameter, truth_value(parameter) * factor)
    }
  }
  paste(row, collapse = ",")
}, character(1))
candidate_preamble <- template[seq_len(header_index - 1L)]
candidate_preamble <- sub(
  "^# num_samples = .*$", paste0("# num_samples = ", nrow(grid)),
  candidate_preamble
)
candidate_csv <- file.path(output_path, "likelihood-candidates.csv")
writeLines(c(candidate_preamble, template[header_index], candidate_rows), candidate_csv)
write.csv(grid, file.path(output_path, "likelihood-candidate-manifest.csv"),
          row.names = FALSE)

sim_model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core_simulate.stan"),
  include_paths = stan_path, cpp_options = list(stan_threads = TRUE)
)
compact_model <- cmdstan_model(
  file.path(stan_path, "takeup_struct_main_core_compact_gq.stan"),
  include_paths = stan_path, cpp_options = list(stan_threads = TRUE)
)

# Evaluate the true reader-facing objects at exactly 1.25 km.
truth_gq_data <- data
truth_gq_data$roc_distances[[6L]] <- truth_gq_data$roc_distances[[6L]] * 2.5
truth_gq <- compact_model$generate_quantities(
  fitted_params = truth_csv, data = truth_gq_data, seed = seed,
  output_dir = output_path, output_basename = "truth-estimands",
  threads_per_chain = 1
)
truth_draw <- as_draws_df(truth_gq$draws(variables = c(
  "core_compact_sm_rescaled", "core_compact_signal_lambda"
)))
truth_values <- data.frame(
  object = c("delta", "lambda", "sigma_u", "psi",
             paste0("multiplier_", c("control", "ink", "calendar", "bracelet"))),
  variable = c("dist_beta_v[1]", "base_mu_rep", "raw_u_sd[1]",
               "wtp_value_utility", paste0("core_compact_sm_rescaled[1,", 1:4, "]")),
  truth = c(
    draws$`dist_beta_v[1]`[medoid_index], draws$base_mu_rep[medoid_index],
    draws$`raw_u_sd[1]`[medoid_index], draws$wtp_value_utility[medoid_index],
    -as.numeric(unlist(truth_draw[1, paste0(
      "core_compact_sm_rescaled[1,", 1:4, "]"
    )], use.names = FALSE))
  )
)
multiplier_truth <- truth_values$truth[startsWith(truth_values$object, "multiplier_")]
truth_values <- bind_rows(
  truth_values,
  data.frame(
    object = "contrast_no_signal_minus_any_signal",
    variable = "derived",
    truth = (multiplier_truth[1] + multiplier_truth[3] -
             multiplier_truth[2] - multiplier_truth[4]) / 2
  )
)
write.csv(truth_values, file.path(output_path, "truth-values.csv"), row.names = FALSE)

sim_data <- data
sim_data$core_sim_zero_lambda <- 0L
sim_data$core_sim_lambda_log_deviation <- rep(0, data$num_treatments)
manifest <- vector("list", replicates)
for (replicate in seq_len(replicates)) {
  replicate_seed <- seed + replicate
  replicate_path <- file.path(dataset_path, sprintf("rep-%03d", replicate))
  dir.create(replicate_path, recursive = TRUE, showWarnings = FALSE)
  gq <- sim_model$generate_quantities(
    fitted_params = truth_csv, data = sim_data, seed = replicate_seed,
    output_dir = replicate_path, output_basename = "simulated-outcomes",
    threads_per_chain = 1
  )
  sim_draw <- as_draws_df(gq$draws(variables = c(
    "core_sim_takeup", "core_sim_num_knows_1ord",
    "core_sim_num_knows_2ord", "core_sim_wtp_response",
    "core_sim_gift_choice"
  )))
  extract_integer <- function(prefix, size) {
    as.integer(unlist(sim_draw[1, paste0(prefix, "[", seq_len(size), "]")],
                      use.names = FALSE))
  }
  refit_data <- data
  refit_data$takeup <- extract_integer("core_sim_takeup", data$num_obs)
  refit_data$num_knows_1ord <- extract_integer(
    "core_sim_num_knows_1ord", data$num_beliefs_obs
  )
  refit_data$num_knows_2ord <- extract_integer(
    "core_sim_num_knows_2ord", data$num_beliefs_obs
  )
  refit_data$wtp_response <- extract_integer(
    "core_sim_wtp_response", data$num_wtp_obs
  )
  refit_data$gift_choice <- extract_integer(
    "core_sim_gift_choice", data$num_wtp_obs
  )
  if (length(refit_data$optim_distances) == 1L &&
      is.null(dim(refit_data$optim_distances))) {
    refit_data$optim_distances <- array(refit_data$optim_distances, dim = 1L)
  }
  data_json <- file.path(replicate_path, "fit-data.json")
  write_stan_json(refit_data, data_json)
  manifest[[replicate]] <- data.frame(
    task_id = replicate, replicate = replicate, seed = replicate_seed,
    data_json = normalizePath(data_json), status = "pending",
    stringsAsFactors = FALSE
  )
}
manifest <- bind_rows(manifest)
write.csv(manifest, file.path(output_path, "benchmark-recovery-manifest.csv"),
          row.names = FALSE)

source_files <- c(workspace, fit_csvs, truth_csv, candidate_csv)
provenance <- data.frame(
  file = normalizePath(source_files), md5 = unname(tools::md5sum(source_files)),
  git_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[1],
  distance_definition = "assigned", replicates = replicates, seed = seed,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write.csv(provenance, file.path(output_path, "provenance.csv"), row.names = FALSE)
message("Prepared ", replicates, " benchmark recovery datasets in ", output_path)
