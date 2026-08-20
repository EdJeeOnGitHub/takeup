#!/usr/bin/env Rscript

root <- "archive/code"
files <- list.files(root, recursive = TRUE, full.names = FALSE)
files <- files[!files %in% c("README.md", "manifest.csv")]
original_path <- sub("^(root-r|root-shell|policy-v1)/", "", files)
original_path <- sub("^scratch/", "scratch/", original_path)
category <- sub("/.*$", "", files)
reason <- ifelse(
  category == "policy-v1",
  "Legacy dense policy workflow; current workflow is scripts/policy/.",
  ifelse(
    category == "scratch",
    "Historical or exploratory scratch file outside the supported paper and appendix workflows.",
    "Historical root entrypoint outside the supported production workflow."
  )
)
replacement <- ifelse(
  grepl("postprocess", basename(files)),
  "scripts/structural/ and scripts/structural/render-paper.R",
  ""
)
manifest <- data.frame(
  original_path = original_path,
  archived_path = file.path(root, files),
  category = category,
  reason = reason,
  replacement = replacement,
  last_pre_archive_commit = "9f05899",
  archived_utc = "2026-08-20",
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$original_path), ]
write.csv(manifest, file.path(root, "manifest.csv"), row.names = FALSE)
cat("Recorded", nrow(manifest), "archived files.\n")
