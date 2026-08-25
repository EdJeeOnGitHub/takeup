#!/usr/bin/env Rscript

# Run memory-efficient generated quantities against existing fit-105 draws.
# This does not sample or alter posterior draws.

args <- commandArgs(trailingOnly = TRUE)

option_value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}
has_flag <- function(name) name %in% args

model_key <- option_value("--model")
input_path <- option_value(
  "--input-path",
  "/project/akaring/takeup-data/data/stan_analysis_data"
)
output_path <- option_value(
  "--output-path",
  "/project/akaring/takeup-data/temp-data/compact-sm-gq"
)
stan_path <- option_value("--stan-path", "stan_models")
parallel_chains <- as.integer(option_value("--parallel-chains", "4"))
threads_per_chain <- as.integer(option_value("--threads-per-chain", "1"))
schema_debug <- has_flag("--schema-debug")
sm_evaluation_distance_m <- as.numeric(
  option_value("--sm-evaluation-distance-m", "500")
)
fit_csv_option <- option_value("--fit-csvs")
chains <- as.integer(strsplit(option_value("--chains", "1,2,3,4"), ",", fixed = TRUE)[[1L]])
if (
  length(chains) == 0L ||
    anyNA(chains) ||
    any(!chains %in% 1:4) ||
    anyDuplicated(chains)
) {
  stop("--chains must be a comma-separated subset of 1,2,3,4.", call. = FALSE)
}
parallel_chains <- min(parallel_chains, length(chains))

model_manifest <- list(
  indiv_dist_community_fp_indiv_vis = list(
    model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS",
    workspace_fit = 104L,
    source_model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP",
    stan_file = "takeup_struct_indiv_vis_sm_compact.stan",
    stan_model_name = "takeup_struct_private_info_model"
  ),
  indiv_dist_indiv_fp = list(
    model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP",
    workspace_fit = 105L,
    source_model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP",
    stan_file = "takeup_struct_indiv_fp_sm_compact.stan",
    stan_model_name = "takeup_struct_no_generated_quantities_model"
  ),
  no_outliers = list(
    model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
    workspace_fit = 1250L,
    source_model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
    stan_file = "takeup_struct_indiv_fp_sm_compact.stan",
    stan_model_name = "takeup_struct_model"
  ),
  correct_observability = list(
    model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS",
    workspace_fit = 1251L,
    source_model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS",
    stan_file = "takeup_struct_indiv_fp_sm_compact.stan",
    stan_model_name = "takeup_struct_model"
  ),
  second_order_observability = list(
    model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB",
    workspace_fit = 1252L,
    source_model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB",
    stan_file = "takeup_struct_indiv_fp_sm_compact.stan",
    stan_model_name = "takeup_struct_model"
  )
)

if (is.null(model_key) || !model_key %in% names(model_manifest)) {
  stop(
    "--model must be one of: ",
    paste(names(model_manifest), collapse = ", "),
    call. = FALSE
  )
}
spec <- model_manifest[[model_key]]

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(purrr)
  library(rlang)
})

workspace_path <- file.path(
  input_path,
  sprintf("dist_fit%d.RData", spec$workspace_fit)
)
if (!file.exists(workspace_path)) {
  stop("Saved fit workspace not found: ", workspace_path, call. = FALSE)
}

fit_env <- new.env(parent = globalenv())
load(workspace_path, envir = fit_env)
if (!all(c("models", "stan_data") %in% ls(fit_env))) {
  stop("Workspace lacks models or stan_data: ", workspace_path, call. = FALSE)
}
if (!spec$source_model %in% names(fit_env$models)) {
  stop(
    "Workspace ", workspace_path, " lacks model configuration ",
    spec$source_model, ".", call. = FALSE
  )
}

model_info <- fit_env$models[[spec$source_model]]
stan_data_preprocess <- model_info$stan_data_preprocess %||% identity
model_info$stan_data_preprocess <- NULL
model_info$model_file <- spec$stan_file

gq_data <- stan_data_preprocess(fit_env$stan_data) |>
  map_at(
    c("cluster_treatment_map", "beliefs_ate_pairs"),
    \(x) mutate(x, across(.cols = everything(), .fns = as.integer)) |>
      as.matrix()
  ) |>
  list_modify(!!!model_info) |>
  map_if(is.factor, as.integer)

# Both fit-105 individual robustness runs used --num-mix-groups=1. The fit-104
# workspace supplies the otherwise-identical community-level Stan data but was
# saved with two mixture groups.
gq_data$num_dist_group_mix <- 1L

# Backward-compatible defaults for belief-data fields added after fit 105.
# These fits used observed belief rows only and did not use row-level beliefs
# to replace the cluster-level reputational return.
gq_data$use_belief_row_cluster_mu_rep <-
  gq_data$use_belief_row_cluster_mu_rep %||% 0L
gq_data$num_belief_rows_by_cluster <-
  gq_data$num_belief_rows_by_cluster %||% rep.int(0L, gq_data$num_clusters)
gq_data$belief_observed <-
  gq_data$belief_observed %||% rep.int(1L, gq_data$num_beliefs_obs)

# Match fit_model() in R/structural/legacy-utils.R, while dropping run controls that
# are not Stan data.
gq_data$control <- NULL
gq_data$analysis_data <- NULL
gq_data <- discard(
  gq_data,
  \(x) is.function(x) || is.character(x) || is.null(x)
)

if (length(gq_data$roc_distances) < 26L) {
  stop("Stan data does not contain the 2500m ROC grid point.", call. = FALSE)
}
if (!is.finite(sm_evaluation_distance_m) || sm_evaluation_distance_m < 0) {
  stop("--sm-evaluation-distance-m must be a finite nonnegative distance.")
}
gq_data$roc_distances[[6L]] <-
  gq_data$roc_distances[[6L]] * sm_evaluation_distance_m / 500

fit_csv <- if (is.null(fit_csv_option)) {
  file.path(
    input_path,
    sprintf("dist_fit105_%s-%d.csv", spec$model, chains)
  )
} else {
  explicit_csv <- strsplit(fit_csv_option, ",", fixed = TRUE)[[1L]]
  if (length(explicit_csv) != length(chains)) {
    stop(
      "--fit-csvs must provide one CSV for every requested chain.",
      call. = FALSE
    )
  }
  explicit_csv
}
if (!all(file.exists(fit_csv))) {
  stop(
    "Missing fitted-parameter CSV(s): ",
    paste(basename(fit_csv[!file.exists(fit_csv)]), collapse = ", "),
    call. = FALSE
  )
}

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
stan_file <- file.path(stan_path, spec$stan_file)
message("Compiling compact GQ model: ", stan_file)
gq_model <- cmdstan_model(
  stan_file,
  include_paths = stan_path,
  cpp_options = list(stan_threads = TRUE),
  stanc_options = list("name" = spec$stan_model_name),
  force_recompile = TRUE
)

if (schema_debug) {
  schema_fit <- gq_model$sample(
    data = gq_data,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 0,
    iter_sampling = 1,
    fixed_param = TRUE,
    init = 0,
    output_dir = output_path,
    output_basename = paste0("schema_", spec$model),
    threads_per_chain = threads_per_chain,
    refresh = 1
  )
  parameter_names <- names(gq_model$variables()$parameters)
  base_name <- function(x) sub("\\[.*$", "", x)
  expected <- schema_fit$metadata()$model_params
  expected <- expected[base_name(expected) %in% parameter_names]
  observed_fit <- as_cmdstan_fit(fit_csv)
  observed <- observed_fit$metadata()$model_params
  observed <- observed[base_name(observed) %in% parameter_names]
  message("Expected parameter columns: ", length(expected))
  message("Observed parameter columns: ", length(observed))
  message(
    "Expected but absent: ",
    paste(setdiff(expected, observed), collapse = ", ")
  )
  message(
    "Observed but unexpected: ",
    paste(setdiff(observed, expected), collapse = ", ")
  )
  quit(status = 0L)
}

message(
  "Running compact GQ for ", spec$model, " with ",
  parallel_chains, " parallel chains and ", threads_per_chain,
  " thread(s) per chain."
)
gq_fit <- gq_model$generate_quantities(
  fitted_params = fit_csv,
  data = gq_data,
  output_dir = output_path,
  output_basename = paste0("compact_sm_fit105_", spec$model),
  parallel_chains = parallel_chains,
  threads_per_chain = threads_per_chain
)

output_files <- gq_fit$output_files()
if (length(output_files) != length(chains) || !all(file.exists(output_files))) {
  stop("Compact GQ did not produce every requested output chain.", call. = FALSE)
}
message("Compact GQ outputs:")
message(paste(output_files, collapse = "\n"))
