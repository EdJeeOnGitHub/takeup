#!/usr/bin/env Rscript

stage <- file.path("build", "replication-package")
if (dir.exists(stage)) unlink(stage, recursive = TRUE)
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
  "R/README.md", "scripts/README.md", "stan_models/README.md",
  "docs/repository-layout.md", "docs/path-migration.md",
  "replication/README.md",
  "replication/data-manifest.csv", "replication/paper-artifact-contract.csv",
  list.files("R", recursive = TRUE, full.names = TRUE),
  unlist(lapply(
    file.path("scripts", c(
      "balance", "checks", "policy", "reduced-form", "shared",
      "structural", "workflow"
    )),
    list.files, recursive = TRUE, full.names = TRUE
  )),
  list.files("tests/benchmarks", recursive = TRUE, full.names = TRUE),
  list.files("tests/smoke", recursive = TRUE, full.names = TRUE),
  list.files("replication/inputs/policy", full.names = TRUE),
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
