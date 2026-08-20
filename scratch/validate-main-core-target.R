#!/usr/bin/env Rscript

# Numerically compare the frozen generic main model and the specialized core
# model at identical unconstrained parameter points. The only intended target
# aggregation replaces repeated Bernoulli terms with the unnormalized
# cluster-binomial kernel, which should be the identical Stan sampling target.

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")
option_value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

input_path <- option_value(
  "--input-path",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-smoke-input"
)
model_name <- option_value(
  "--model",
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS"
)
stan_path <- option_value("--stan-path", "stan_models_fit105")
generic_file <- option_value("--generic-file", "takeup_struct.stan")
core_file <- option_value("--core-file", "takeup_struct_main_core.stan")
output_path <- option_value("--output-path", "temp/main-core-validation")
cmdstan_path_option <- option_value(
  "--cmdstan-path",
  Sys.getenv("CMDSTAN", unset = "/home/edjee/.cmdstan/cmdstan-2.33.1")
)
lp_tolerance <- as.numeric(option_value("--lp-tolerance", "1e-6"))
gradient_tolerance <- as.numeric(
  option_value("--gradient-tolerance", "1e-6")
)

options(cmdstanr_no_ver_check = TRUE)
Sys.setenv(
  CMDSTAN = cmdstan_path_option,
  CMDSTANR_NO_VER_CHECK = "TRUE"
)
suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(purrr)
  library(rlang)
})
set_cmdstan_path(cmdstan_path_option)

workspace_path <- file.path(input_path, "dist_fit105.RData")
fit_env <- new.env(parent = globalenv())
load(workspace_path, envir = fit_env)
if (!model_name %in% names(fit_env$models)) {
  stop("Workspace lacks model configuration: ", model_name, call. = FALSE)
}

model_info <- fit_env$models[[model_name]]
stan_data_preprocess <- model_info$stan_data_preprocess %||% identity
model_info$stan_data_preprocess <- NULL
model_info$model_file <- NULL

sample_data <- stan_data_preprocess(fit_env$stan_data) |>
  map_at(
    c("cluster_treatment_map", "beliefs_ate_pairs"),
    \(x) mutate(x, across(everything(), as.integer)) |> as.matrix()
  ) |>
  list_modify(!!!model_info) |>
  map_if(is.factor, as.integer)
sample_data$num_dist_group_mix <- 1L
sample_data$use_belief_row_cluster_mu_rep <-
  sample_data$use_belief_row_cluster_mu_rep %||% 0L
if (length(sample_data$num_belief_rows_by_cluster) != sample_data$num_clusters) {
  sample_data$num_belief_rows_by_cluster <- tabulate(
    sample_data$obs_cluster_id[sample_data$beliefs_obs_index],
    nbins = sample_data$num_clusters
  )
}
if (length(sample_data$belief_observed) != sample_data$num_beliefs_obs) {
  sample_data$belief_observed <- rep.int(1L, sample_data$num_beliefs_obs)
}
if (sample_data$num_optim_distances == 1L) {
  sample_data$optim_distances <- array(sample_data$optim_distances, dim = 1L)
}
sample_data$core_cluster_weight <- rep(1, sample_data$num_clusters)
sample_data$core_visibility_prior_multiplier <- 1
sample_data$core_wtp_mu_prior_sd <- 2
sample_data$use_core_cluster_shock <- 0L
sample_data$core_cluster_shock_sd_prior <- 0.1
sample_data$core_lambda_structure <- 0L
sample_data$core_lambda_log_ratio_sd_prior <- 0.25
sample_data$core_profile_group_lambda <- 0L
sample_data$core_profile_group_log_ratio <- 0
sample_data$core_gq_override_lambda <- 0L
sample_data$core_gq_lambda_override <- rep(0, sample_data$num_treatments)
type_mixture <- main_core_type_mixture_data(0L, 5, 12)
sample_data$core_type_distribution <- 0L
sample_data$core_student_t_df <- type_mixture$df
sample_data$core_type_scale_sq <- type_mixture$scale_sq
sample_data$core_type_mixture_components <- type_mixture$components
sample_data$core_type_mixture_precision <- type_mixture$precision
sample_data$core_type_mixture_weight <- type_mixture$weight
sample_data$core_observation_model <- 0L
sample_data <- modifyList(sample_data, main_core_empty_peer_response_data())
sample_data$core_peer_link_audit <- NULL
# Legacy workspaces did not retain this mapping. With unit weights its values
# cannot affect the target, so a valid placeholder preserves exact equivalence.
sample_data$wtp_cluster_id <- sample_data$wtp_cluster_id %||%
  rep.int(1L, sample_data$num_wtp_obs)
sample_data$control <- NULL
sample_data$analysis_data <- NULL
sample_data <- discard(
  sample_data,
  \(x) is.function(x) || is.character(x) || is.null(x)
)

# Reproduce takeup_transformed_data.stan's included monitored sample.
all_clusters <- seq_len(sample_data$num_clusters)
included_clusters <- if (sample_data$num_excluded_clusters > 0L) {
  setdiff(all_clusters, sample_data$excluded_clusters)
} else {
  all_clusters
}
included_obs <- which(sample_data$obs_cluster_id %in% included_clusters)
monitored_obs <- which(sample_data$is_name_matched == 0L)
included_monitored_obs <- intersect(monitored_obs, included_obs)
included_cluster_id <- sample_data$obs_cluster_id[included_monitored_obs]
included_takeup <- sample_data$takeup[included_monitored_obs]

cluster_n <- tabulate(included_cluster_id, nbins = sample_data$num_clusters)
cluster_y <- tabulate(
  included_cluster_id[included_takeup == 1L],
  nbins = sample_data$num_clusters
)
binomial_constant <- sum(lchoose(cluster_n[included_clusters], cluster_y[included_clusters]))

# Test the Bernoulli-to-binomial identity directly at several heterogeneous
# cluster-probability vectors. This is the only intentional likelihood rewrite.
set.seed(20260802L)
probability_tests <- list(
  constant = rep(0.37, sample_data$num_clusters),
  deterministic = seq(0.05, 0.95, length.out = sample_data$num_clusters),
  random = runif(sample_data$num_clusters, 0.02, 0.98)
)
likelihood_identity <- imap_dfr(probability_tests, function(probability, label) {
  bernoulli_lp <- sum(dbinom(
    included_takeup,
    size = 1L,
    prob = probability[included_cluster_id],
    log = TRUE
  ))
  binomial_lp <- sum(dbinom(
    cluster_y[included_clusters],
    size = cluster_n[included_clusters],
    prob = probability[included_clusters],
    log = TRUE
  ))
  tibble(
    probability_test = label,
    bernoulli_lp = bernoulli_lp,
    binomial_lp = binomial_lp,
    omitted_binomial_constant = binomial_constant,
    identity_residual =
      (binomial_lp - binomial_constant) - bernoulli_lp
  )
})

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
generic_model <- cmdstan_model(
  file.path(stan_path, generic_file),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
core_model <- cmdstan_model(
  file.path(stan_path, core_file),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)
# Let CmdStanR perform model-aware JSON serialization (notably preserving
# length-one arrays), then reuse that exact data file for high-precision direct
# CmdStan calls below.
data_probe <- generic_model$diagnose(
  data = sample_data,
  init = 0,
  seed = 20260801L,
  error = 1,
  output_dir = output_path
)
data_json <- data_probe$data_file()

diagnose_at <- function(model, init, seed, label) {
  diagnostic_output <- system2(
    model$exe_file(),
    c(
      "diagnose",
      "test=gradient",
      "epsilon=1e-6",
      "error=1",
      "data",
      paste0("file=", data_json),
      paste0("init=", init),
      "random",
      paste0("seed=", seed),
      "output",
      "sig_figs=18"
    ),
    stdout = TRUE,
    stderr = TRUE,
    env = c("STAN_NUM_THREADS=1")
  )
  status <- attr(diagnostic_output, "status") %||% 0L
  if (status != 0L) {
    stop(
      "CmdStan diagnose failed for ", label, ":\n",
      paste(diagnostic_output, collapse = "\n"),
      call. = FALSE
    )
  }

  lp_line <- grep("Log probability=", diagnostic_output, value = TRUE)
  if (length(lp_line) != 1L) {
    stop("Could not parse log probability for ", label, call. = FALSE)
  }
  lp <- as.numeric(sub(".*Log probability=", "", lp_line))

  gradient_lines <- grep("^[[:space:]]+[0-9]+[[:space:]]", diagnostic_output, value = TRUE)
  gradients <- vapply(
    strsplit(trimws(gradient_lines), "[[:space:]]+"),
    \(fields) as.numeric(fields[[3L]]),
    numeric(1)
  )
  if (length(gradients) == 0L || anyNA(gradients)) {
    stop("Could not parse gradients for ", label, call. = FALSE)
  }
  list(lp = lp, gradients = gradients)
}

test_points <- list(
  list(label = "zero", init = 0, seed = 20260802L),
  list(label = "random_small_1", init = 0.01, seed = 20260803L),
  list(label = "random_small_2", init = 0.01, seed = 20260804L)
)

results <- map_dfr(test_points, function(point) {
  point_label <- point$label
  point_seed <- point$seed
  message("Diagnosing target at point: ", point_label)
  generic <- diagnose_at(
    generic_model,
    point$init,
    point_seed,
    paste0("generic/", point_label)
  )
  core <- diagnose_at(
    core_model,
    point$init,
    point_seed,
    paste0("core/", point_label)
  )

  generic_lp <- generic$lp
  core_lp <- core$lp
  generic_gradients <- generic$gradients
  core_gradients <- core$gradients
  if (length(generic_gradients) != length(core_gradients)) {
    stop("Gradient lengths differ at ", point_label, call. = FALSE)
  }
  tibble(
    point = point_label,
    seed = point_seed,
    # CmdStan 2.33's diagnostic summary rounds this line to six significant
    # figures regardless of output sig_figs; retain it only as a coarse check.
    generic_lp_coarse = generic_lp,
    core_lp_coarse = core_lp,
    max_abs_gradient_difference = max(
      abs(core_gradients - generic_gradients)
    )
  )
})

report_path <- file.path(output_path, "main-core-target-equivalence.csv")
write.csv(results, report_path, row.names = FALSE)
identity_report_path <- file.path(
  output_path,
  "main-core-likelihood-identity.csv"
)
write.csv(likelihood_identity, identity_report_path, row.names = FALSE)
print(results, n = Inf)
print(likelihood_identity, n = Inf)

if (any(abs(likelihood_identity$identity_residual) > lp_tolerance)) {
  stop(
    "Bernoulli/binomial likelihood identity failed. See ",
    identity_report_path,
    call. = FALSE
  )
}
if (any(results$max_abs_gradient_difference > gradient_tolerance)) {
  stop("Gradient equivalence test failed. See ", report_path, call. = FALSE)
}
message(
  "Target-equivalence tests passed: ",
  report_path,
  " and ",
  identity_report_path
)
