#!/usr/bin/env Rscript

# Sample a fit-105 slim structural model while saving only parameters.
# Large deterministic arrays remain local to the Stan model.

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")

option_value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

model_name <- option_value(
  "--model",
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP"
)
input_path <- option_value(
  "--input-path",
  "/project/akaring/takeup-data/data/stan_analysis_data"
)
output_path <- option_value(
  "--output-path",
  "/project/akaring/takeup-data/data/stan_analysis_data"
)
data_json <- option_value("--data-json")
stan_path <- option_value("--stan-path", "stan_models")
stan_file_name <- option_value(
  "--stan-file",
  "takeup_struct_indiv_fp_slim.stan"
)
cmdstan_path_option <- option_value(
  "--cmdstan-path",
  Sys.getenv("CMDSTAN_PATH", unset = "")
)
output_basename <- option_value(
  "--output-basename",
  paste0("dist_fit106_", model_name, "_SLIM")
)
chains <- as.integer(option_value("--chains", "4"))
parallel_chains <- as.integer(option_value("--parallel-chains", as.character(chains)))
threads_per_chain <- as.integer(option_value("--threads-per-chain", "2"))
iter_warmup <- as.integer(option_value("--iter-warmup", "1000"))
iter_sampling <- as.integer(option_value("--iter-sampling", "1000"))
adapt_delta <- as.numeric(option_value("--adapt-delta", "0.999"))
max_treedepth <- as.integer(option_value("--max-treedepth", "12"))
metric <- option_value("--metric", "dense_e")
seed <- as.integer(option_value("--seed", "20260724"))
refresh <- as.integer(option_value("--refresh", "25"))
init_files_option <- option_value("--init-files")
fit_model_override <- option_value("--fit-model")
fit_dist_override <- option_value("--fit-dist-model")
fit_beliefs_override <- option_value("--fit-beliefs-model")
fit_wtp_override <- option_value("--fit-wtp-model")
cluster_weight_file <- option_value("--cluster-weight-file")
use_core_cluster_shock <- as.integer(
  option_value("--use-core-cluster-shock", "0")
)
core_cluster_shock_sd_prior <- as.numeric(
  option_value("--core-cluster-shock-sd-prior", "0.1")
)
core_lambda_structure <- as.integer(
  option_value("--core-lambda-structure", "0")
)
core_lambda_log_ratio_sd_prior <- as.numeric(
  option_value("--core-lambda-log-ratio-sd-prior", "0.25")
)
core_profile_group_lambda <- as.integer(
  option_value("--core-profile-group-lambda", "0")
)
core_profile_group_log_ratio <- as.numeric(
  option_value("--core-profile-group-log-ratio", "0")
)
core_type_distribution <- as.integer(
  option_value("--core-type-distribution", "0")
)
core_student_t_df <- as.numeric(option_value("--core-student-t-df", "5"))
core_student_t_components <- as.integer(
  option_value("--core-student-t-components", "12")
)
core_observation_model <- as.integer(
  option_value("--core-observation-model", "0")
)
core_recognition_structure <- as.integer(
  option_value("--core-recognition-structure", "0")
)
core_report_structure <- as.integer(
  option_value("--core-report-structure", "0")
)

stopifnot(
  chains >= 1L,
  parallel_chains >= 1L,
  parallel_chains <= chains,
  threads_per_chain >= 1L,
  iter_warmup >= 0L,
  iter_sampling >= 1L,
  adapt_delta > 0,
  adapt_delta < 1,
  max_treedepth >= 1L,
  metric %in% c("diag_e", "dense_e")
)
if (!use_core_cluster_shock %in% 0:1) {
  stop("--use-core-cluster-shock must be 0 or 1.", call. = FALSE)
}
if (!is.finite(core_cluster_shock_sd_prior) ||
    core_cluster_shock_sd_prior <= 0) {
  stop("--core-cluster-shock-sd-prior must be positive.", call. = FALSE)
}
if (!core_lambda_structure %in% 0:2) {
  stop("--core-lambda-structure must be 0, 1, or 2.", call. = FALSE)
}
if (!is.finite(core_lambda_log_ratio_sd_prior) ||
    core_lambda_log_ratio_sd_prior <= 0) {
  stop("--core-lambda-log-ratio-sd-prior must be positive.", call. = FALSE)
}
if (!core_profile_group_lambda %in% 0:1 ||
    (core_profile_group_lambda == 1L && core_lambda_structure != 1L)) {
  stop("Grouped-lambda profiling requires structure 1.", call. = FALSE)
}
if (!core_type_distribution %in% 0:1 || !is.finite(core_student_t_df) ||
    core_student_t_df <= 2 || core_student_t_components < 2L) {
  stop("Student-t controls require distribution 0/1, df > 2, and components >= 2.",
       call. = FALSE)
}
if (!core_observation_model %in% 0:2) {
  stop("--core-observation-model must be 0, 1, or 2.", call. = FALSE)
}
if (!core_recognition_structure %in% 0:2 ||
    !core_report_structure %in% 0:2) {
  stop("Observation structures must be 0, 1, or 2.", call. = FALSE)
}
if (core_observation_model == 0L &&
    (core_recognition_structure != 0L || core_report_structure != 0L)) {
  stop("Observation restrictions require an asymmetric model.", call. = FALSE)
}
if (core_observation_model == 2L && core_recognition_structure == 2L) {
  stop("The unconditional model cannot condition recognition out.",
       call. = FALSE)
}

init_value <- if (is.null(init_files_option)) {
  0
} else {
  init_files <- strsplit(init_files_option, ",", fixed = TRUE)[[1L]]
  if (length(init_files) != chains || !all(file.exists(init_files))) {
    stop(
      "--init-files must provide one existing JSON file per chain.",
      call. = FALSE
    )
  }
  # Mode JSONs can contain exact lower-bound values.  Nudge those values into
  # the parameter space before passing them to HMC.
  sanitized_init_files <- file.path(
    output_path,
    sprintf(
      "%s-sanitized-init-%02d.json",
      output_basename,
      seq_along(init_files)
    )
  )
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  for (index in seq_along(init_files)) {
    init <- jsonlite::read_json(init_files[index], simplifyVector = TRUE)
    init <- main_core_nudge_init_boundaries(init)
    # jsonlite simplifies a one-element JSON array to a scalar. Restore the
    # singleton arrays used by the core Stan parameter schema.
    for (parameter in intersect(
      names(init), c(
        "raw_u_sd", "dist_beta_v", "core_cluster_shock_sd",
        "core_lambda_group_log_ratio_raw"
      )
    )) {
      if (is.null(dim(init[[parameter]])) && length(init[[parameter]]) == 1L) {
        init[[parameter]] <- array(init[[parameter]], dim = 1L)
      }
    }
    cmdstanr::write_stan_json(init, sanitized_init_files[index])
  }
  sanitized_init_files
}

# Compute nodes do not need GitHub access to check for newer CmdStan releases.
options(cmdstanr_no_ver_check = TRUE)
Sys.setenv(CMDSTANR_NO_VER_CHECK = "TRUE")
old_timeout <- getOption("timeout")
options(timeout = 3)
if (nzchar(cmdstan_path_option)) {
  # CmdStanR discovers this environment variable during namespace loading.
  Sys.setenv(CMDSTAN = cmdstan_path_option)
}
message("Loading CmdStanR (remote version check disabled).")
suppressPackageStartupMessages(library(cmdstanr))
options(timeout = old_timeout)
message("Loading data-preparation packages.")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(rlang)
})
if (nzchar(cmdstan_path_option)) {
  if (!dir.exists(cmdstan_path_option)) {
    stop("CmdStan path not found: ", cmdstan_path_option, call. = FALSE)
  }
  set_cmdstan_path(cmdstan_path_option)
}

workspace_path <- file.path(input_path, "dist_fit105.RData")
if (!file.exists(workspace_path)) {
  stop("Saved fit workspace not found: ", workspace_path, call. = FALSE)
}

fit_env <- new.env(parent = globalenv())
load(workspace_path, envir = fit_env)
if (!all(c("models", "stan_data") %in% ls(fit_env))) {
  stop("Workspace lacks models or stan_data: ", workspace_path, call. = FALSE)
}
if (!model_name %in% names(fit_env$models)) {
  stop("Workspace lacks model configuration: ", model_name, call. = FALSE)
}

model_info <- fit_env$models[[model_name]]
stan_data_preprocess <- model_info$stan_data_preprocess %||% identity
model_info$stan_data_preprocess <- NULL
model_info$model_file <- stan_file_name

sample_data <- stan_data_preprocess(fit_env$stan_data) |>
  map_at(
    c("cluster_treatment_map", "beliefs_ate_pairs"),
    \(x) mutate(x, across(.cols = everything(), .fns = as.integer)) |>
      as.matrix()
  ) |>
  list_modify(!!!model_info) |>
  map_if(is.factor, as.integer)

# Match the original fit-105 invocation.
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

# Minimal main-model extensions. Older saved workspaces do not retain the WTP
# cluster mapping; that is harmless for an unweighted fit but must be rejected
# for cluster-weighted refits.
sample_data$core_cluster_weight <- rep(1, sample_data$num_clusters)
sample_data$use_core_cluster_shock <- use_core_cluster_shock
sample_data$core_cluster_shock_sd_prior <- core_cluster_shock_sd_prior
sample_data$core_lambda_structure <- core_lambda_structure
sample_data$core_lambda_log_ratio_sd_prior <-
  core_lambda_log_ratio_sd_prior
sample_data$core_profile_group_lambda <- core_profile_group_lambda
sample_data$core_profile_group_log_ratio <- core_profile_group_log_ratio
sample_data$core_gq_override_lambda <- 0L
sample_data$core_gq_lambda_override <- rep(0, sample_data$num_treatments)
type_mixture <- main_core_student_t_mixture(
  core_student_t_df, core_student_t_components
)
sample_data$core_type_distribution <- core_type_distribution
sample_data$core_student_t_df <- type_mixture$df
sample_data$core_type_scale_sq <- type_mixture$scale_sq
sample_data$core_type_mixture_components <- type_mixture$components
sample_data$core_type_mixture_precision <- type_mixture$precision
sample_data$core_type_mixture_weight <- type_mixture$weight
sample_data$core_observation_model <- core_observation_model
sample_data$core_recognition_structure <- core_recognition_structure
sample_data$core_report_structure <- core_report_structure
peer_data <- if (core_observation_model == 0L) {
  main_core_empty_peer_response_data()
} else {
  main_core_prepare_peer_response_data(
    sample_data,
    project_root = ".",
    write_audit = file.path(
      output_path,
      paste0(output_basename, "-peer-link-audit.csv")
    )
  )
}
peer_data$core_peer_link_audit <- NULL
sample_data <- modifyList(sample_data, peer_data)
if (is.null(sample_data$wtp_cluster_id)) {
  if (!is.null(cluster_weight_file)) {
    stop(
      "The workspace lacks wtp_cluster_id; regenerate the main workspace ",
      "before running weighted refits.",
      call. = FALSE
    )
  }
  warning(
    "Workspace lacks wtp_cluster_id; assigning WTP rows to cluster 1. ",
    "This is target-equivalent only because all cluster weights are one."
  )
  sample_data$wtp_cluster_id <- rep.int(1L, sample_data$num_wtp_obs)
}
if (!is.null(cluster_weight_file)) {
  if (!file.exists(cluster_weight_file)) {
    stop("Cluster-weight file not found: ", cluster_weight_file, call. = FALSE)
  }
  weight_data <- read.csv(cluster_weight_file, stringsAsFactors = FALSE)
  if (!identical(names(weight_data), c("cluster_id", "weight"))) {
    stop(
      "Cluster-weight CSV must have exactly: cluster_id,weight.",
      call. = FALSE
    )
  }
  if (nrow(weight_data) != sample_data$num_clusters ||
      anyDuplicated(weight_data$cluster_id) ||
      !setequal(weight_data$cluster_id, seq_len(sample_data$num_clusters)) ||
      any(!is.finite(weight_data$weight)) ||
      any(weight_data$weight < 0)) {
    stop("Invalid cluster-weight CSV.", call. = FALSE)
  }
  weight_data <- weight_data[order(weight_data$cluster_id), ]
  sample_data$core_cluster_weight <- weight_data$weight
}
sample_data$control <- NULL
sample_data$analysis_data <- NULL
sample_data <- discard(
  sample_data,
  \(x) is.function(x) || is.character(x) || is.null(x)
)
if (!is.null(data_json)) {
  if (core_observation_model > 0L) {
    stop("Asymmetric observability requires rebuilding data, not --data-json.",
         call. = FALSE)
  }
  if (!file.exists(data_json)) stop("Data JSON not found: ", data_json, call. = FALSE)
  patched_data_json <- tempfile(fileext = ".json")
  on.exit(unlink(patched_data_json), add = TRUE)
  sample_data <- main_core_patch_stan_json_scalars(
    data_json,
    patched_data_json,
    c(
      core_lambda_structure = core_lambda_structure,
      core_lambda_log_ratio_sd_prior = core_lambda_log_ratio_sd_prior,
      core_profile_group_lambda = core_profile_group_lambda,
      core_profile_group_log_ratio = core_profile_group_log_ratio,
      core_gq_override_lambda = 0L,
      core_type_distribution = core_type_distribution,
      core_student_t_df = core_student_t_df,
      core_type_scale_sq = type_mixture$scale_sq,
      core_type_mixture_components = core_student_t_components,
      core_observation_model = core_observation_model,
      core_recognition_structure = core_recognition_structure,
      core_report_structure = core_report_structure
    )
  )
}

parse_binary_override <- function(value, option_name) {
  if (is.null(value)) return(NULL)
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || !parsed %in% 0:1) {
    stop(option_name, " must be 0 or 1.", call. = FALSE)
  }
  parsed
}

likelihood_overrides <- list(
  fit_model_to_data = parse_binary_override(
    fit_model_override,
    "--fit-model"
  ),
  fit_dist_model_to_data = parse_binary_override(
    fit_dist_override,
    "--fit-dist-model"
  ),
  fit_beliefs_model_to_data = parse_binary_override(
    fit_beliefs_override,
    "--fit-beliefs-model"
  ),
  fit_wtp_model_to_data = parse_binary_override(
    fit_wtp_override,
    "--fit-wtp-model"
  )
)
for (override_name in names(likelihood_overrides)) {
  override_value <- likelihood_overrides[[override_name]]
  if (!is.null(override_value)) {
    sample_data[[override_name]] <- override_value
  }
}

expected_num_clusters <- option_value("--expected-num-clusters")
if (!is.null(expected_num_clusters) &&
    sample_data$num_clusters != as.integer(expected_num_clusters)) {
  stop(
    "Expected ", expected_num_clusters, " clusters; found ",
    sample_data$num_clusters,
    ".",
    call. = FALSE
  )
}

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
stan_file <- file.path(stan_path, stan_file_name)
message("Compiling slim structural model: ", stan_file)
model <- cmdstan_model(
  stan_file,
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE)
)

message(
  "Sampling ", chains, " chains (", parallel_chains, " parallel), ",
  iter_warmup, " warmup + ", iter_sampling, " retained draws, ",
  threads_per_chain, " thread(s)/chain, metric=", metric,
  ", adapt_delta=", adapt_delta, "."
)

fit <- model$sample(
  data = sample_data,
  chains = chains,
  parallel_chains = parallel_chains,
  threads_per_chain = threads_per_chain,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  save_warmup = FALSE,
  thin = 1,
  # The default avoids a historical init closure whose global enums were not
  # serialized. Production runs can instead supply posterior warm-start JSON.
  init = init_value,
  seed = seed,
  metric = metric,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  output_dir = output_path,
  output_basename = output_basename,
  refresh = refresh
)

message("Output files:")
message(paste(fit$output_files(), collapse = "\n"))
print(fit$diagnostic_summary())
