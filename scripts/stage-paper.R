#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- sub("^--spec=", "", args[grepl("^--spec=", args)])
if (!length(value)) value <- Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned")
strict <- "--strict" %in% args
source("R/distance-spec.R")
specification <- takeup_distance_spec(value[[1L]])

system2("Rscript", c("--no-save", "--no-restore",
                     "scripts/build-paper-artifact-registry.R"))
registry <- read.csv("build/manifest/paper-artifacts.csv",
                     stringsAsFactors = FALSE)
contract <- read.csv("replication/paper-artifact-contract.csv",
                     stringsAsFactors = FALSE, na.strings = character())
registry <- merge(registry, contract, by = "artifact", all.x = TRUE,
                  sort = FALSE, suffixes = c("", "_contract"))
stage <- file.path("build", "paper", specification)
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
file.copy("ref-reports/ECM ReStud.tex", file.path(stage, "ECM ReStud.tex"),
          overwrite = TRUE)

unresolved <- character()
for (index in seq_len(nrow(registry))) {
  destination <- file.path(stage, registry$artifact[index])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  generated_matches <- list.files(
    file.path("build", "work", specification), recursive = TRUE,
    full.names = TRUE,
    pattern = paste0("^", gsub("([.(){}+*?^$|\\[\\]\\\\])", "\\\\\\1",
                                basename(registry$artifact[index])), "$")
  )
  generated <- file.path("build", specification, registry$source_path[index])
  source <- if (registry$default_contract[index] == "generated" &&
                length(generated_matches)) {
    generated_matches[[1L]]
  } else if (registry$default_contract[index] == "generated" &&
             file.exists(generated)) {
    generated
  } else if (registry$default_contract[index] == "frozen" &&
             file.exists(registry$deposit_path[index])) {
    registry$deposit_path[index]
  } else {
    registry$source_path[index]
  }
  if (!length(source) || is.na(source) || !file.exists(source)) {
    unresolved <- c(unresolved, registry$artifact[index])
  } else if (!file.copy(source, destination, overwrite = TRUE)) {
    unresolved <- c(unresolved, registry$artifact[index])
  }
}
writeLines(unresolved, file.path(stage, "missing-artifacts.txt"))
if (length(unresolved) && strict) {
  stop(length(unresolved), " contracted paper artifacts could not be staged:\n",
       paste(unresolved, collapse = "\n"))
} else if (length(unresolved)) {
  warning(length(unresolved), " contracted paper artifacts could not be staged.")
} else if (nzchar(Sys.which("latexmk"))) {
  old_dir <- setwd(stage)
  on.exit(setwd(old_dir), add = TRUE)
  status <- system2("latexmk", c("-pdf", "-interaction=nonstopmode",
                                  shQuote("ECM ReStud.tex")), stdout = TRUE,
                    stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    pdf <- "ECM ReStud.pdf"
    if (file.exists(pdf)) {
      warning("LaTeX reported errors, but produced ", file.path(stage, pdf),
              ". The current manuscript source has pre-existing TeX errors; ",
              "inspect its log before release.")
    } else {
      stop("Paper compilation failed in ", stage, call. = FALSE)
    }
  }
}
cat(stage, "\n")
