#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) default else substring(hit[[1L]], nchar(prefix) + 1L)
}

specification <- value("--spec")
output_dir <- value("--output-dir")
basename <- value(
  "--basename", "dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)
warmup <- as.integer(value("--iter-warmup", "1000"))
sampling <- as.integer(value("--iter-sampling", "1000"))
seed <- as.integer(value("--seed", "20260820"))
workspace <- value(
  "--workspace",
  Sys.getenv(
    "TAKEUP_STRUCTURAL_WORKSPACE",
    "build/structural-workspace/main-core-input.RData"
  )
)
if (is.null(specification) || is.null(output_dir)) quit(status = 1L)
if (!file.exists(workspace)) quit(status = 1L)
workspace_sha256 <- strsplit(
  system2("sha256sum", workspace, stdout = TRUE),
  "[[:space:]]+"
)[[1L]][[1L]]

chains <- file.path(output_dir, sprintf("%s-%d.csv", basename, 1:4))
manifest_path <- file.path(output_dir, paste0(basename, "-manifest.csv"))
diagnostic_path <- file.path(output_dir, paste0(basename, "-diagnostics.csv"))
required <- c(chains, manifest_path, diagnostic_path)
info <- file.info(required)
if (any(is.na(info$size)) || any(info$size == 0L)) quit(status = 1L)

manifest <- tryCatch(read.csv(manifest_path, stringsAsFactors = FALSE),
                     error = function(e) NULL)
diagnostics <- tryCatch(read.csv(diagnostic_path, stringsAsFactors = FALSE),
                        error = function(e) NULL)
if (is.null(manifest) || nrow(manifest) != 1L ||
    is.null(diagnostics) || nrow(diagnostics) != 4L) quit(status = 1L)

valid_manifest <- identical(manifest$distance_definition[[1L]], specification) &&
  manifest$chains[[1L]] == 4L &&
  manifest$iter_warmup[[1L]] == warmup &&
  manifest$iter_sampling[[1L]] == sampling &&
  manifest$seed[[1L]] == seed &&
  "workspace_sha256" %in% names(manifest) &&
  identical(manifest$workspace_sha256[[1L]], workspace_sha256)
valid_diagnostics <-
  all(diagnostics$num_divergent == 0L) &&
  all(diagnostics$num_max_treedepth == 0L) &&
  all(is.finite(diagnostics$ebfmi)) && all(diagnostics$ebfmi > 0.3)

if (!valid_manifest || !valid_diagnostics) quit(status = 1L)
cat("Validated cached four-chain structural fit in", output_dir, "\n")
