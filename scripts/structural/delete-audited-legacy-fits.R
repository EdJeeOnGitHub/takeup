#!/usr/bin/env Rscript

# Delete only the exact legacy CmdStan CSVs approved by the streamlined-fit
# equivalence audit.  A SHA-256 receipt is written before any file is removed.

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit)) substring(hit, nchar(prefix) + 1L) else default
}

analysis_root <- normalizePath(value(
  "--analysis-root",
  "/project/akaring/takeup-data/data/stan_analysis_data"
), mustWork = TRUE)
audit_path <- value(
  "--audit-path",
  file.path(
    analysis_root, "streamlined-active-robustness", "audit",
    "legacy-deletion-manifest.csv"
  )
)
receipt_path <- value(
  "--receipt-path",
  file.path(dirname(audit_path), "legacy-deletion-receipt.csv")
)
confirm <- identical(value("--confirm-delete"), "YES")
if (!confirm) {
  stop("Refusing deletion without --confirm-delete=YES.", call. = FALSE)
}
if (!file.exists(audit_path)) stop("Missing audit manifest: ", audit_path)

manifest <- read.csv(audit_path, stringsAsFactors = FALSE)
required <- c(
  "spec_id", "path", "exists", "bytes", "replacement_root",
  "delete_approved_by_audit"
)
missing <- setdiff(required, names(manifest))
if (length(missing)) stop("Audit manifest is missing: ", paste(missing, collapse = ", "))
if (nrow(manifest) != 16L || anyDuplicated(manifest$path)) {
  stop("Expected exactly 16 unique audited legacy paths.", call. = FALSE)
}
if (!all(manifest$exists) || !all(manifest$delete_approved_by_audit)) {
  stop("Every manifest row must exist and be approved by the audit.", call. = FALSE)
}

resolved_parent <- vapply(
  manifest$path, function(path) normalizePath(dirname(path), mustWork = TRUE),
  character(1)
)
if (!all(resolved_parent == analysis_root)) {
  stop("Every deletion target must be directly beneath the analysis root.")
}
allowed_pattern <- paste0(
  "^dist_fit10[56]_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_",
  "(INDIV_DIST_COMMUNITY_FP_INDIV_VIS|INDIV_DIST_INDIV_FP|NO_OUTLIERS|SOB)-",
  "[1-4][.]csv$"
)
if (!all(grepl(allowed_pattern, basename(manifest$path)))) {
  stop("At least one target is outside the four audited fit families.")
}
actual_bytes <- file.info(manifest$path)$size
if (any(is.na(actual_bytes)) || !identical(as.numeric(actual_bytes), as.numeric(manifest$bytes))) {
  stop("Target sizes no longer match the audited manifest.")
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path)
  sub("[[:space:]].*$", "", output[[1L]])
}

receipt <- transform(
  manifest,
  sha256 = vapply(manifest$path, sha256_file, character(1)),
  deletion_started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  deleted = FALSE,
  exists_after = TRUE
)
write.csv(receipt, receipt_path, row.names = FALSE, quote = TRUE)

deleted <- file.remove(manifest$path)
receipt$deleted <- deleted
receipt$exists_after <- file.exists(manifest$path)
receipt$deletion_finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
write.csv(receipt, receipt_path, row.names = FALSE, quote = TRUE)

if (!all(deleted) || any(receipt$exists_after)) {
  stop("At least one audited legacy file could not be removed; inspect ", receipt_path)
}
message(
  "Deleted ", nrow(receipt), " audited legacy CSVs (",
  sum(receipt$bytes), " bytes). Receipt: ", receipt_path
)
