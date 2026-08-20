#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- startsWith(args, prefix)
  if (any(hit)) return(sub(prefix, "", args[which(hit)[1L]], fixed = TRUE))
  index <- match(name, args)
  if (!is.na(index) && index < length(args)) return(args[index + 1L])
  default
}

owner <- option_value("--owner")
generated_root <- option_value("--generated-root")
output <- option_value("--output")
contract_path <- option_value(
  "--contract", "replication/paper-artifact-contract.csv"
)
if (is.null(owner) || is.null(generated_root) || is.null(output)) {
  stop("--owner, --generated-root, and --output are required.", call. = FALSE)
}
if (!dir.exists(generated_root)) {
  stop("Generated root does not exist: ", generated_root, call. = FALSE)
}

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}
contract <- read.csv(contract_path, stringsAsFactors = FALSE)
contract <- contract[contract$workflow_owner == owner, , drop = FALSE]
generated <- list.files(generated_root, recursive = TRUE, full.names = TRUE)
generated <- generated[file.info(generated)$isdir %in% FALSE]
comparison_name <- function(path) {
  sub("-1\\.pdf$", ".pdf", basename(path))
}

rows <- lapply(seq_len(nrow(contract)), function(index) {
  artifact <- contract$artifact[index]
  matches <- generated[comparison_name(generated) == comparison_name(artifact)]
  generated_path <- if (length(matches) == 1L) matches else NA_character_
  approved_path <- contract$deposit_path[index]
  approved_hash <- sha256(approved_path)
  generated_hash <- sha256(generated_path)
  data.frame(
    artifact = artifact,
    generated_path = generated_path,
    approved_path = approved_path,
    status = if (is.na(generated_path)) {
      "not_rendered"
    } else if (identical(generated_hash, approved_hash)) {
      "byte_identical"
    } else {
      "changed_review_required"
    },
    generated_md5 = generated_hash,
    approved_md5 = approved_hash,
    stringsAsFactors = FALSE
  )
})
comparison <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, output, row.names = FALSE)
message(
  "Compared ", sum(comparison$status != "not_rendered"), " rendered artifacts; ",
  sum(comparison$status == "changed_review_required"), " require review."
)
