#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (!length(hit)) default else substring(hit[[1L]], nchar(name) + 2L)
}
specification <- value("--spec", "assigned")
stage <- value("--stage", file.path("build", "paper-candidate", specification))
required <- file.path(stage, c(
  "ECM ReStud.pdf", "structural-robustness-appendix.pdf",
  "artifact-comparison.csv", "build-manifest.csv", "numeric-text-audit.csv",
  "NUMERIC-TODOS.md", "candidate-status.json", "missing-artifacts.txt"
))
missing_files <- required[!file.exists(required)]
if (length(missing_files)) stop("Missing candidate outputs:\n", paste(missing_files, collapse = "\n"))
if (length(readLines(file.path(stage, "missing-artifacts.txt"), warn = FALSE))) {
  stop("Candidate contains unresolved artifacts.")
}
manifest <- read.csv(file.path(stage, "build-manifest.csv"), stringsAsFactors = FALSE)
if (manifest$distance_definition[[1L]] != specification ||
    manifest$missing_artifacts[[1L]] != 0L ||
    !manifest$manuscript_compiled[[1L]] || !manifest$appendix_compiled[[1L]]) {
  stop("Candidate build manifest failed validation.")
}
artifacts <- read.csv(file.path(stage, "artifact-comparison.csv"), stringsAsFactors = FALSE)
bad <- artifacts$distance_sensitive & artifacts$selected_contract != "candidate"
if (any(bad)) stop("Distance-sensitive artifacts lack candidate provenance:\n",
                   paste(artifacts$artifact[bad], collapse = "\n"))

rf_marker <- file.path("build", specification, "reduced-form.complete")
if (!file.exists(rf_marker)) stop("Missing reduced-form production marker.")
rf <- read.csv(rf_marker, stringsAsFactors = FALSE)
if (!all(c("distance_definition", "bootstrap_draws") %in% names(rf)) ||
    rf$distance_definition[[1L]] != specification || rf$bootstrap_draws[[1L]] < 500L) {
  stop("Reduced-form output lacks assigned 500-draw provenance.")
}
balance_sections <- c("main", "orig", "fit-ri", "attrition", "monitored-attrition", "sms")
for (section in balance_sections) {
  marker <- file.path("build", specification, paste0("balance-", section, ".complete"))
  if (!file.exists(marker)) stop("Missing balance marker: ", section)
  info <- read.csv(marker, stringsAsFactors = FALSE)
  if (!all(c("distance_definition", "ri_draws", "required_output_sha256") %in% names(info)) ||
      info$distance_definition[[1L]] != specification || info$ri_draws[[1L]] < 500L) {
    stop("Balance section lacks assigned 500-draw provenance: ", section)
  }
}

fit_manifest <- file.path(
  "build", "structural-fit", specification,
  "dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-manifest.csv"
)
if (!file.exists(fit_manifest)) stop("Missing assigned slim structural-fit manifest.")
fit_info <- read.csv(fit_manifest, stringsAsFactors = FALSE)
if (fit_info$distance_definition[[1L]] != specification || fit_info$chains[[1L]] != 4L ||
    fit_info$iter_sampling[[1L]] < 1000L) stop("Structural fit is not a production assigned run.")
diagnostic_path <- sub("-manifest[.]csv$", "-diagnostics.csv", fit_manifest)
if (!file.exists(diagnostic_path)) stop("Missing structural diagnostic audit.")
diagnostics <- read.csv(diagnostic_path, stringsAsFactors = FALSE)
for (field in intersect(c("num_divergent", "num_max_treedepth"), names(diagnostics))) {
  if (any(diagnostics[[field]] > 0L)) stop("Structural diagnostic failure: ", field)
}

hpc_root <- file.path("build", "candidate-hpc", specification, "imported")
for (name in c("candidate-hpc-manifest.csv", "artifact-manifest.csv")) {
  if (!file.exists(file.path(hpc_root, name))) stop("Missing validated HPC import: ", name)
}
hpc_artifacts <- read.csv(file.path(hpc_root, "artifact-manifest.csv"), stringsAsFactors = FALSE)
if (any(hpc_artifacts$distance_definition != specification)) stop("HPC import mixes distance definitions.")
policy_manifest <- file.path(
  hpc_root, "artifacts", "policy", "policy-bootstrap-manifest.csv"
)
if (!file.exists(policy_manifest)) stop("Missing assigned policy parameter manifest.")
policy <- read.csv(policy_manifest, stringsAsFactors = FALSE)
if (policy$distance_definition[[1L]] != specification || policy$replications[[1L]] < 999L) {
  stop("Policy parameters are not the assigned 999-replicate production set.")
}

numeric_audit <- read.csv(file.path(stage, "numeric-text-audit.csv"), stringsAsFactors = FALSE)
if (!all(c("claim_id", "old_value", "status", "context") %in% names(numeric_audit)) ||
    anyDuplicated(numeric_audit$claim_id)) stop("Numerical prose TODO audit is invalid.")
cat("Candidate build is complete; numerical prose remains review-pending.\n")
