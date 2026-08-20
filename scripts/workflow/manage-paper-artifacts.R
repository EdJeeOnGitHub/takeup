#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
initialize <- "--initialize" %in% args
contract_path <- "replication/paper-artifact-contract.csv"
deposit_root <- "replication/paper-artifacts"

sha256 <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path)
  sub("[[:space:]].*$", "", output[[1L]])
}

if (initialize) {
  status <- system2("Rscript", c("--no-save", "--no-restore",
    "scripts/checks/audit-paper-pipeline-coverage.R"))
  if (!identical(status, 0L)) stop("Could not audit manuscript artifacts.")
  x <- read.csv("build/manifest/paper-pipeline-coverage.csv",
                stringsAsFactors = FALSE)

  recovered <- c(
    "misc-figures/comp-dist-plot3-fit86-util-identity-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.pdf" =
      "/home/agent/projects/overleaf/overleaf-takeup-slides/Stanford/comp-dist-plot3-fit86-util-identity-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.pdf",
    "misc-figures/all-clusters-map-1.pdf" =
      "/home/agent/projects/overleaf/overleaf-takeup-slides/all-clusters-map-1.pdf"
  )
  reconstructed <- c("timeline-tikz.tex", "anchor-diagram.tex")
  if ("default_contract" %in% names(x)) {
    freeze <- x$default_contract == "frozen"
  } else {
    freeze <- !(x$pipeline_coverage %in%
      c("covered_make_target", "covered_versioned_source"))
    x$default_contract <- ifelse(
      x$pipeline_coverage == "covered_make_target", "generated",
      ifelse(x$pipeline_coverage == "covered_versioned_source", "static", "frozen")
    )
  }
  x$deposit_path <- ifelse(freeze, file.path(deposit_root, x$artifact), "")

  for (i in which(freeze)) {
    destination <- x$deposit_path[[i]]
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (x$artifact[[i]] %in% reconstructed && file.exists(destination)) next
    source <- x$source_path[[i]]
    if (is.na(source) || !file.exists(source)) {
      source <- unname(recovered[x$artifact[[i]]])
    }
    if (!length(source) || is.na(source) || !file.exists(source)) {
      stop("No approved artifact available for ", x$artifact[[i]])
    }
    if (!file.exists(destination) && !file.copy(source, destination)) {
      stop("Could not deposit ", x$artifact[[i]])
    }
  }

  full_target <- c(
    balance = "balance",
    reduced_form = "reduced-form",
    structural_postprocess = "structural-postprocess",
    optimal_policy = "optimal-policy policy-tables",
    design_simulation = "design-paper-full",
    versioned_source_asset = "paper-assets",
    manual_or_static = "paper-assets"
  )
  assigned_target <- unname(full_target[x$workflow_owner])
  assigned_target[freeze & x$workflow_owner %in% c("reduced_form", "balance")] <-
    "auxiliary-paper"
  contract <- data.frame(
    artifact = x$artifact,
    workflow_owner = x$workflow_owner,
    default_contract = x$default_contract,
    source_path = ifelse(is.na(x$source_path), "", x$source_path),
    deposit_path = x$deposit_path,
    sha256 = ifelse(freeze, vapply(x$deposit_path, sha256, character(1)), ""),
    full_make_target = assigned_target,
    provenance_note = ifelse(
      x$artifact %in% names(recovered), "Recovered from the presentation repository.",
      ifelse(x$artifact %in% reconstructed,
        "Reconstructed from the manuscript design description; review visually.",
        ifelse(freeze, "Frozen current paper-approved artifact.",
          ifelse(x$default_contract == "generated",
            "Regenerated in the isolated Make/targets build.",
            "Version-controlled source asset.")))
    ),
    stringsAsFactors = FALSE
  )
  write.csv(contract, contract_path, row.names = FALSE, na = "")
  cat("Initialized", sum(freeze), "frozen artifacts and wrote", contract_path, "\n")
}

if (!file.exists(contract_path)) {
  stop("Artifact contract missing. Run this script once with --initialize.")
}
contract <- read.csv(contract_path, stringsAsFactors = FALSE,
                     na.strings = character())
frozen <- contract$default_contract == "frozen"
missing <- contract$deposit_path[frozen & !file.exists(contract$deposit_path)]
if (length(missing)) stop("Missing frozen artifacts:\n", paste(missing, collapse = "\n"))
actual <- vapply(contract$deposit_path[frozen], sha256, character(1))
bad <- which(actual != contract$sha256[frozen])
if (length(bad)) {
  stop("Frozen artifact checksum mismatch:\n",
       paste(contract$artifact[frozen][bad], collapse = "\n"))
}
static <- contract$source_path[contract$default_contract == "static"]
missing_static <- static[!file.exists(static)]
if (length(missing_static)) {
  stop("Missing versioned source assets:\n", paste(missing_static, collapse = "\n"))
}
cat("Validated", sum(frozen), "frozen and", length(static), "static paper assets.\n")
