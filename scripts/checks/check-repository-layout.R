#!/usr/bin/env Rscript

fail <- function(...) stop(..., call. = FALSE)

is_checkout <- file.exists(".git")
required <- c(
  "R/README.md", "scripts/README.md", "stan_models/README.md",
  "docs/repository-layout.md", "docs/path-migration.md"
)
if (is_checkout) required <- c(
  required, "hpc/README.md", "scratch/README.md",
  "archive/code/README.md", "archive/code/manifest.csv"
)
missing <- required[!file.exists(required)]
if (length(missing)) fail("Required layout files missing: ", paste(missing, collapse = ", "))

if (is_checkout) {
  tracked <- system2("git", c("ls-files"), stdout = TRUE)
  root_entrypoints <- setdiff(
    grep("^[^/]+\\.(R|sh)$", tracked, value = TRUE),
    "_targets.R"
  )
  if (length(root_entrypoints)) {
    fail("R or shell entry points remain at repository root: ", paste(root_entrypoints, collapse = ", "))
  }
}

if (is_checkout) {
  archive_files <- list.files("archive/code", recursive = TRUE, full.names = TRUE)
  archive_files <- archive_files[!basename(archive_files) %in% c("README.md", "manifest.csv")]
  manifest <- read.csv("archive/code/manifest.csv", stringsAsFactors = FALSE)
  unrecorded <- setdiff(archive_files, manifest$archived_path)
  stale <- setdiff(manifest$archived_path, archive_files)
  if (length(unrecorded) || length(stale)) {
    fail(
      "Archive manifest mismatch. Unrecorded: ", paste(unrecorded, collapse = ", "),
      "; missing: ", paste(stale, collapse = ", ")
    )
  }
}

production_files <- c(
  "Makefile", "_targets.R",
  list.files("R", recursive = TRUE, full.names = TRUE),
  list.files("scripts", recursive = TRUE, full.names = TRUE),
  list.files("hpc", recursive = TRUE, full.names = TRUE)
)
production_files <- production_files[file.info(production_files)$isdir %in% FALSE]
production_files <- setdiff(production_files, c(
  "scripts/checks/check-repository-layout.R",
  "scripts/workflow/build-code-archive-manifest.R"
))
scratch_refs <- unlist(lapply(production_files, function(path) {
  lines <- readLines(path, warn = FALSE)
  hits <- grep("(^|[[:space:]\\\"'])scratch/", lines, perl = TRUE)
  if (!length(hits)) return(character())
  paste0(path, ":", hits, ":", trimws(lines[hits]))
}), use.names = FALSE)
scratch_refs <- scratch_refs[!grepl("archive/code/scratch/", scratch_refs, fixed = TRUE)]
if (length(scratch_refs)) {
  fail("Production files refer to repository scratch/:\n", paste(scratch_refs, collapse = "\n"))
}

cat("Repository layout checks passed.\n")
