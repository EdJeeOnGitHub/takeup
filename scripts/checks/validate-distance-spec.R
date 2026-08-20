#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- sub("^--spec=", "", args[grepl("^--spec=", args)])
if (!length(value)) value <- Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned")

source("R/distance/spec.R")
specification <- takeup_distance_spec(value[[1L]])
crosswalk <- takeup_distance_crosswalk()
output_dir <- file.path("build", specification, "audit")
paths <- takeup_write_distance_audit(crosswalk, output_dir, specification)

analysis_crosswalk <- dplyr::filter(crosswalk, .data$in_main_analysis)
cat("Validated", nrow(analysis_crosswalk), "clusters using", specification, "Close/Far.\n")
cat("Changed classifications:", sum(analysis_crosswalk$switched), "\n")
cat(paste(paths, collapse = "\n"), "\n")
