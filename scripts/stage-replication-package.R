#!/usr/bin/env Rscript

stage <- file.path("build", "replication-package")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)

collect_stan_dependencies <- function(entrypoints, stan_dir = "stan_models") {
  queue <- entrypoints
  found <- character()
  while (length(queue)) {
    path <- queue[[1L]]
    queue <- queue[-1L]
    if (path %in% found || !file.exists(path)) next
    found <- c(found, path)
    lines <- readLines(path, warn = FALSE)
    includes <- sub(
      ".*#include[[:space:]]+([^[:space:]]+).*", "\\1",
      grep("#include[[:space:]]+", lines, value = TRUE)
    )
    queue <- c(queue, file.path(stan_dir, basename(includes)))
  }
  unique(found)
}

stan_sources <- collect_stan_dependencies(c(
  "stan_models/takeup_struct_main_core.stan",
  "stan_models/takeup_struct_main_core_compact_gq.stan"
))

artifact_contract <- read.csv(
  "replication/paper-artifact-contract.csv",
  stringsAsFactors = FALSE,
  na.strings = character()
)
contract_files <- c(
  artifact_contract$deposit_path[
    artifact_contract$default_contract == "frozen"
  ],
  artifact_contract$source_path[
    artifact_contract$default_contract == "static"
  ]
)

include <- c(
  "README.md", "Makefile", "_targets.R", "renv.lock", "takeup.Rproj",
  "R/distance-spec.R", "R/pipeline.R", "replication/README.md",
  "replication/data-manifest.csv", "replication/paper-artifact-contract.csv",
  "scripts/validate-distance-spec.R",
  "scripts/prepare-structural-distance-data.R",
  "scripts/run-structural-fit.sh", "scripts/build-paper-artifact-registry.R",
  "scripts/audit-paper-pipeline-coverage.R",
  "scripts/build-stan-artifact-inventory.R",
  "scripts/check-replication-build.R", "scripts/stage-paper.R",
  "scripts/manage-paper-artifacts.R",
  "scripts/stage-replication-package.R",
  "scripts/test-distance-spec.R",
  "scratch/reduced-form-bootstrap.R", "scratch/reduced-form-setup.R",
  "scratch/reduced-form-functions.R", "scratch/bench-rf-direct-wls.R",
  "scratch/main-core-data.R",
  "scratch/sample-slim-individual-fp.R",
  "scratch/generate-main-core-compact-gq.R",
  "ref-reports/ECM ReStud.tex",
  stan_sources, contract_files
)
include <- unique(include)
missing <- include[!file.exists(include)]
if (length(missing)) stop("Package inputs missing: ", paste(missing, collapse = ", "))
for (path in include) {
  destination <- file.path(stage, path)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  file.copy(path, destination, overwrite = TRUE)
}

files <- list.files(stage, recursive = TRUE, full.names = TRUE)
files <- files[file.info(files)$isdir %in% FALSE]
checksums <- tools::md5sum(files)
manifest <- data.frame(
  path = sub(paste0("^", normalizePath(stage), "/"), "", normalizePath(files)),
  bytes = file.info(files)$size,
  md5 = unname(checksums)
)
write.csv(manifest, file.path(stage, "checksums.csv"), row.names = FALSE)
cat(stage, "\n")
