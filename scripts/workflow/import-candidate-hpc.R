#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (!length(hit)) default else substring(hit[[1L]], nchar(name) + 2L)
}
specification <- value("--spec", "assigned")
source("R/distance/spec.R")
specification <- takeup_distance_spec(specification)
export_root <- value("--export-root")
bundle <- value("--bundle")
output_root <- value("--output-root")
if (any(!nzchar(c(export_root, bundle, output_root)))) stop("Export, bundle, and output roots are required.")

expected <- read.csv(file.path(export_root, "candidate-hpc-manifest.csv"), stringsAsFactors = FALSE)
received <- read.csv(file.path(bundle, "candidate-hpc-manifest.csv"), stringsAsFactors = FALSE)
if (!identical(expected, received)) stop("HPC bundle manifest does not match the exported request.")
if (!identical(received$distance_definition[[1L]], specification)) stop("Wrong distance definition in HPC bundle.")
artifacts <- read.csv(file.path(bundle, "artifact-manifest.csv"), stringsAsFactors = FALSE)
required_columns <- c("path", "sha256", "workflow_id", "distance_definition")
if (!all(required_columns %in% names(artifacts)) || !nrow(artifacts)) stop("Invalid or empty artifact manifest.")
if (any(artifacts$distance_definition != specification) || anyDuplicated(artifacts$path) ||
    any(grepl("(^|/)\\.\\.(/|$)|^/", artifacts$path))) stop("Invalid artifact provenance or path.")
jobs <- read.csv(file.path(export_root, "job-manifest.csv"), stringsAsFactors = FALSE)
missing_workflows <- setdiff(jobs$workflow_id[jobs$required], unique(artifacts$workflow_id))
if (length(missing_workflows)) stop("HPC bundle lacks required workflows: ",
                                    paste(missing_workflows, collapse = ", "))
expected_paths <- unique(unlist(strsplit(
  paste(jobs$expected_artifacts[jobs$required], collapse = ";"), ";", fixed = TRUE
)))
expected_paths <- expected_paths[nzchar(expected_paths)]
missing_expected <- setdiff(expected_paths, artifacts$path)
if (length(missing_expected)) stop("HPC bundle lacks expected artifacts:\n",
                                   paste(missing_expected, collapse = "\n"))
paths <- file.path(bundle, "artifacts", artifacts$path)
if (any(!file.exists(paths))) stop("HPC bundle is missing declared artifacts.")
actual <- sub("[[:space:]].*$", "", vapply(paths, function(path) system2("sha256sum", path, stdout = TRUE)[[1L]], character(1)))
if (any(actual != artifacts$sha256)) stop("HPC artifact checksum mismatch.")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
existing <- if (dir.exists(file.path(output_root, "artifacts"))) {
  list.files(file.path(output_root, "artifacts"), recursive = TRUE)
} else character()
extra <- setdiff(existing, artifacts$path)
if (length(extra)) stop("Import destination contains stale undeclared artifacts; use a clean candidate build root.")
for (index in seq_along(paths)) {
  destination <- file.path(output_root, "artifacts", artifacts$path[[index]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(paths[[index]], destination, overwrite = TRUE)) stop("Could not import ", artifacts$path[[index]])
}
file.copy(file.path(bundle, "artifact-manifest.csv"), output_root, overwrite = TRUE)
file.copy(file.path(bundle, "candidate-hpc-manifest.csv"), output_root, overwrite = TRUE)
cat(output_root, "\n")
