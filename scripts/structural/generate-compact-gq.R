#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(purrr)
  library(rlang)
})

workspace_path <- option_value(
  "--workspace", "data/stan_analysis_data/dist_fit104.RData"
)
data_json <- option_value("--data-json")
model_name <- option_value(
  "--model", "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)
stan_path <- option_value("--stan-path", "stan_models")
stan_file <- option_value(
  "--stan-file", "takeup_struct_main_core_compact_gq.stan"
)
fit_csvs <- strsplit(
  option_value("--fit-csvs", ""), ",", fixed = TRUE
)[[1L]]
output_path <- option_value("--output-path", "temp/main-core-compact-gq")
output_basename <- option_value("--output-basename")
use_cluster_shock <- as.integer(
  option_value("--use-core-cluster-shock", "0")
)
shock_prior <- as.numeric(
  option_value("--core-cluster-shock-sd-prior", "0.1")
)
lambda_structure <- as.integer(
  option_value("--core-lambda-structure", "0")
)
lambda_prior <- as.numeric(
  option_value("--core-lambda-log-ratio-sd-prior", "0.25")
)
profile_group_lambda <- as.integer(
  option_value("--core-profile-group-lambda", "0")
)
profile_group_log_ratio <- as.numeric(
  option_value("--core-profile-group-log-ratio", "0")
)
type_distribution <- as.integer(option_value("--core-type-distribution", "0"))
student_t_df <- as.numeric(option_value("--core-student-t-df", "5"))
student_t_components <- as.integer(
  option_value("--core-student-t-components", "12")
)
observation_model <- as.integer(
  option_value("--core-observation-model", "0")
)
recognition_structure <- as.integer(
  option_value("--core-recognition-structure", "0")
)
report_structure <- as.integer(
  option_value("--core-report-structure", "0")
)
report_arm_dist_hierarchical <- as.integer(
  option_value("--core-report-arm-dist-hierarchical", "0")
)
report_arm_dist_prior_scale <- as.numeric(
  option_value("--core-report-arm-dist-prior-scale", "0.25")
)
visibility_prior_multiplier <- as.numeric(
  option_value("--core-visibility-prior-multiplier", "1")
)
wtp_mu_prior_sd <- as.numeric(option_value("--core-wtp-mu-prior-sd", "2"))
optional_numeric <- function(name) {
  value <- option_value(name)
  if (is.null(value)) NULL else as.numeric(value)
}
prior_overrides <- list(
  mu_rep_sd = optional_numeric("--mu-rep-sd"),
  dist_beta_v_sd = optional_numeric("--dist-beta-v-sd"),
  raw_u_sd_alpha = optional_numeric("--raw-u-sd-alpha"),
  raw_u_sd_beta = optional_numeric("--raw-u-sd-beta"),
  beta_intercept_sd = optional_numeric("--beta-intercept-sd"),
  beta_ink_effect_sd = optional_numeric("--beta-ink-effect-sd"),
  beta_calendar_effect_sd = optional_numeric("--beta-calendar-effect-sd"),
  beta_bracelet_effect_sd = optional_numeric("--beta-bracelet-effect-sd"),
  wtp_value_utility_mean = optional_numeric("--wtp-value-utility-mean"),
  wtp_value_utility_sd = optional_numeric("--wtp-value-utility-sd"),
  lnorm_wtp_value_utility_prior = optional_numeric(
    "--lnorm-wtp-value-utility-prior"
  )
)
peer_audit_path <- option_value("--peer-audit-path")
cluster_weight_file <- option_value("--cluster-weight-file")
distance_definition <- option_value(
  "--distance-definition", Sys.getenv("TAKEUP_DISTANCE_SPEC", "realized")
)
distance_definition <- takeup_distance_spec(distance_definition)
if (!is.null(data_json) && distance_definition != "realized") {
  stop(
    "Assigned Close/Far requires rebuilding structural data from --workspace; ",
    "an existing --data-json cannot be relabeled safely.", call. = FALSE
  )
}
threads <- as.integer(option_value("--threads", "2"))
parallel_chains <- as.integer(option_value("--parallel-chains", "4"))
cmdstan_path_option <- option_value(
  "--cmdstan-path", Sys.getenv("CMDSTAN_PATH", unset = "")
)
force_recompile <- as.integer(option_value("--force-recompile", "0")) == 1L
sm_evaluation_distance_m <- optional_numeric("--sm-evaluation-distance-m")

if (length(fit_csvs) < 1L || any(!nzchar(fit_csvs)) ||
    any(!file.exists(fit_csvs))) {
  stop("--fit-csvs must list existing fitted-parameter CSVs.", call. = FALSE)
}
if (nzchar(cmdstan_path_option)) {
  set_cmdstan_path(cmdstan_path_option)
}
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE")
options(cmdstanr_no_ver_check = TRUE)

data <- if (!is.null(data_json)) {
  if (observation_model > 0L) {
    stop("Asymmetric observability requires rebuilding data, not --data-json.",
         call. = FALSE)
  }
  if (!file.exists(data_json)) stop("Data JSON not found: ", data_json, call. = FALSE)
  patched_data_json <- tempfile(fileext = ".json")
  on.exit(unlink(patched_data_json), add = TRUE)
  main_core_patch_stan_json_scalars(
    data_json,
    patched_data_json,
    c(
      core_lambda_structure = lambda_structure,
      core_lambda_log_ratio_sd_prior = lambda_prior,
      core_profile_group_lambda = profile_group_lambda,
      core_profile_group_log_ratio = profile_group_log_ratio,
      core_gq_override_lambda = 0L,
      core_visibility_prior_multiplier = visibility_prior_multiplier,
      core_wtp_mu_prior_sd = wtp_mu_prior_sd,
      unlist(
        prior_overrides[!vapply(prior_overrides, is.null, logical(1))],
        use.names = TRUE
      )
    )
  )
} else {
  prepare_main_core_data(
    workspace_path,
    model_name,
    weight_file = cluster_weight_file,
    use_cluster_shock = use_cluster_shock,
    cluster_shock_sd_prior = shock_prior,
    lambda_structure = lambda_structure,
    lambda_log_ratio_sd_prior = lambda_prior,
    profile_group_lambda = profile_group_lambda,
    profile_group_log_ratio = profile_group_log_ratio,
    type_distribution = type_distribution,
    student_t_df = student_t_df,
    student_t_components = student_t_components,
    observation_model = observation_model,
    recognition_structure = recognition_structure,
    report_structure = report_structure,
    report_arm_dist_hierarchical = report_arm_dist_hierarchical,
    report_arm_dist_prior_scale = report_arm_dist_prior_scale,
    visibility_prior_multiplier = visibility_prior_multiplier,
    wtp_mu_prior_sd = wtp_mu_prior_sd,
    mu_rep_sd = prior_overrides$mu_rep_sd,
    dist_beta_v_sd = prior_overrides$dist_beta_v_sd,
    raw_u_sd_alpha = prior_overrides$raw_u_sd_alpha,
    raw_u_sd_beta = prior_overrides$raw_u_sd_beta,
    beta_intercept_sd = prior_overrides$beta_intercept_sd,
    beta_ink_effect_sd = prior_overrides$beta_ink_effect_sd,
    beta_calendar_effect_sd = prior_overrides$beta_calendar_effect_sd,
    beta_bracelet_effect_sd = prior_overrides$beta_bracelet_effect_sd,
    wtp_value_utility_mean = prior_overrides$wtp_value_utility_mean,
    wtp_value_utility_sd = prior_overrides$wtp_value_utility_sd,
    lnorm_wtp_value_utility_prior =
      prior_overrides$lnorm_wtp_value_utility_prior,
    distance_definition = distance_definition,
    peer_audit_path = peer_audit_path
  )
}
if (!is.null(sm_evaluation_distance_m)) {
  if (!is.finite(sm_evaluation_distance_m) || sm_evaluation_distance_m < 0) {
    stop("--sm-evaluation-distance-m must be a finite nonnegative distance.")
  }
  reference_index <- 6L
  reference_distance_m <- 500
  if (length(data$roc_distances) < reference_index ||
      !is.finite(data$roc_distances[[reference_index]]) ||
      data$roc_distances[[reference_index]] <= 0) {
    stop("Stan data lacks the reference 500m ROC grid point.")
  }
  data$roc_distances[[reference_index]] <-
    data$roc_distances[[reference_index]] *
    sm_evaluation_distance_m / reference_distance_m
}
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
model <- cmdstan_model(
  file.path(stan_path, stan_file),
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE),
  force_recompile = force_recompile
)
gq_fit <- model$generate_quantities(
  fitted_params = fit_csvs,
  data = data,
  output_dir = output_path,
  output_basename = if (!is.null(output_basename)) {
    output_basename
  } else if (type_distribution == 1L) {
    paste0("main-core-student-t-df", format(student_t_df, trim = TRUE), "-compact")
  } else if (type_distribution == 2L) {
    "main-core-finite-mixture-compact"
  } else if (observation_model > 0L) {
    paste0(
      "main-core-observation-", observation_model, "-rec-",
      recognition_structure, "-report-", report_structure, "-compact"
    )
  } else if (lambda_structure > 0) {
    paste0(
      "main-core-lambda-", lambda_structure, "-",
      format(lambda_prior, trim = TRUE),
      "-compact"
    )
  } else if (use_cluster_shock) {
    paste0("main-core-shock-", format(shock_prior, trim = TRUE), "-compact")
  } else {
    "main-core-baseline-compact"
  },
  parallel_chains = min(parallel_chains, length(fit_csvs)),
  threads_per_chain = threads
)
gq_files <- gq_fit$output_files()
if (length(gq_files) != length(fit_csvs) ||
    !all(file.exists(gq_files)) || any(file.info(gq_files)$size == 0L)) {
  stop(
    "Generated quantities did not produce one nonempty output per fit CSV.",
    call. = FALSE
  )
}
message(paste(gq_files, collapse = "\n"))
