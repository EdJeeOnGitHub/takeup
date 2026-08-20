#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- sub("^--spec=", "", args[grepl("^--spec=", args)])
if (!length(value)) value <- Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(rlang)
})
source("scratch/main-core-data.R")

specification <- takeup_distance_spec(value[[1L]])
model_name <- Sys.getenv(
  "TAKEUP_STRUCTURAL_MODEL",
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)
workspace <- Sys.getenv("TAKEUP_STRUCTURAL_WORKSPACE", "")
if (!nzchar(workspace)) {
  candidates <- c(
    "data/stan_analysis_data/dist_fit105.RData",
    "data/stan_analysis_data/dist_fit104.RData"
  )
  workspace <- candidates[vapply(candidates, function(path) {
    if (!file.exists(path)) return(FALSE)
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    model_name %in% names(env$models)
  }, logical(1))][1L]
}
if (is.na(workspace) || !file.exists(workspace)) {
  stop("No structural workspace contains model ", model_name, call. = FALSE)
}
fit_env <- new.env(parent = globalenv())
load(workspace, envir = fit_env)
model_info <- fit_env$models[[model_name]]
preprocess <- model_info$stan_data_preprocess %||% identity
model_info$stan_data_preprocess <- NULL
model_info$model_file <- NULL
sample_data <- preprocess(fit_env$stan_data) |>
  map_at(
    c("cluster_treatment_map", "beliefs_ate_pairs"),
    \(x) mutate(x, across(everything(), as.integer)) |> as.matrix()
  ) |>
  list_modify(!!!model_info) |>
  map_if(is.factor, as.integer)
sample_data <- main_core_apply_distance_definition(
  sample_data, specification, project_root = "."
)

analysis_clusters <- sample_data$analysis_data |>
  distinct(cluster_id, cluster.id, assigned.treatment,
           assigned_dist_group, realized_dist_group, analysis_dist_group) |>
  arrange(cluster_id)
analysis_clusters$stan_distance_group_id <-
  sample_data$cluster_assigned_dist_group
analysis_clusters$stan_treatment_group_id <-
  sample_data$cluster_assigned_dist_group_treatment

output_dir <- file.path("build", specification, "structural")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
path <- file.path(output_dir, "structural-distance-input.csv")
write.csv(analysis_clusters, path, row.names = FALSE)
cat(path, "\n")
