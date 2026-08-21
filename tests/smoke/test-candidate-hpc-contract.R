#!/usr/bin/env Rscript

root <- normalizePath(".")
test_root <- file.path(tempdir(), paste0("candidate-hpc-contract-", Sys.getpid()))
export_root <- file.path(test_root, "export")
bundle <- file.path(test_root, "bundle")
imported <- file.path(test_root, "imported")
status <- system2(file.path(R.home("bin"), "Rscript"), c(
  "scripts/workflow/export-candidate-hpc.R", "--spec=assigned",
  paste0("--output-root=", export_root), "--allow-dirty"
))
stopifnot(status == 0L)
dir.create(file.path(bundle, "artifacts"), recursive = TRUE)
file.copy(file.path(export_root, "candidate-hpc-manifest.csv"), bundle)
jobs <- read.csv(file.path(export_root, "job-manifest.csv"), stringsAsFactors = FALSE)
paths <- unique(c(
  file.path("contract-test", paste0(jobs$workflow_id, ".txt")),
  unlist(strsplit(paste(jobs$expected_artifacts, collapse = ";"), ";", fixed = TRUE))
))
paths <- paths[nzchar(paths)]
workflow_for_path <- vapply(paths, function(path) {
  hit <- which(vapply(jobs$expected_artifacts, function(x) {
    path %in% strsplit(x, ";", fixed = TRUE)[[1L]]
  }, logical(1)))
  if (length(hit)) jobs$workflow_id[[hit[[1L]]]] else {
    sub("[.]txt$", "", basename(path))
  }
}, character(1))
for (index in seq_along(paths)) {
  destination <- file.path(bundle, "artifacts", paths[[index]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  writeLines(workflow_for_path[[index]], destination)
}
sha <- vapply(file.path(bundle, "artifacts", paths), function(path) {
  sub("[[:space:]].*$", "", system2("sha256sum", path, stdout = TRUE)[[1L]])
}, character(1))
write.csv(data.frame(
  path = paths, sha256 = sha, workflow_id = workflow_for_path,
  distance_definition = "assigned"
), file.path(bundle, "artifact-manifest.csv"), row.names = FALSE)
status <- system2(file.path(R.home("bin"), "Rscript"), c(
  "scripts/workflow/import-candidate-hpc.R", "--spec=assigned",
  paste0("--export-root=", export_root), paste0("--bundle=", bundle),
  paste0("--output-root=", imported)
))
stopifnot(status == 0L)
stopifnot(all(file.exists(file.path(imported, "artifacts", paths))))
cat("candidate HPC contract smoke test passed\n")
