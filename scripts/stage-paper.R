#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- sub("^--spec=", "", args[grepl("^--spec=", args)])
if (!length(value)) value <- Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned")
source("R/distance-spec.R")
specification <- takeup_distance_spec(value[[1L]])

system2("Rscript", c("--no-save", "--no-restore",
                     "scripts/build-paper-artifact-registry.R"))
registry <- read.csv("build/manifest/paper-artifacts.csv",
                     stringsAsFactors = FALSE)
stage <- file.path("build", "paper", specification)
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
file.copy("ref-reports/ECM ReStud.tex", file.path(stage, "ECM ReStud.tex"),
          overwrite = TRUE)

for (index in seq_len(nrow(registry))) {
  if (registry$status[index] != "resolved") next
  destination <- file.path(stage, registry$artifact[index])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  generated_matches <- list.files(
    file.path("build", "work", specification), recursive = TRUE,
    full.names = TRUE,
    pattern = paste0("^", gsub("([.(){}+*?^$|\\[\\]\\\\])", "\\\\\\1",
                                basename(registry$artifact[index])), "$")
  )
  generated <- file.path("build", specification, registry$source_path[index])
  source <- if (length(generated_matches)) {
    generated_matches[[1L]]
  } else if (file.exists(generated)) {
    generated
  } else {
    registry$source_path[index]
  }
  file.copy(source, destination, overwrite = TRUE)
}
missing <- registry$artifact[registry$status == "missing"]
writeLines(missing, file.path(stage, "missing-artifacts.txt"))
if (length(missing)) {
  warning(length(missing), " paper artifacts are missing from the repository snapshot.")
} else if (nzchar(Sys.which("latexmk"))) {
  old_dir <- setwd(stage)
  on.exit(setwd(old_dir), add = TRUE)
  status <- system2("latexmk", c("-pdf", "-interaction=nonstopmode",
                                  "ECM ReStud.tex"), stdout = TRUE,
                    stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    stop("Paper compilation failed in ", stage, call. = FALSE)
  }
}
cat(stage, "\n")
