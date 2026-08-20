#!/usr/bin/env Rscript

# Classify every active manuscript dependency by the workflow that owns it and
# whether the current Make/targets graph actually regenerates it. Existing
# files are not treated as reproduced merely because stage-paper.R can copy them.

status <- system2(
  "Rscript",
  c("--no-save", "--no-restore", "scripts/build-paper-artifact-registry.R")
)
if (!identical(status, 0L)) stop("Could not build paper artifact registry.")

registry <- read.csv(
  "build/manifest/paper-artifacts.csv",
  stringsAsFactors = FALSE
)
specification <- Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned")

escape_regex <- function(value) {
  gsub("([.(){}+*?^$|\\[\\]\\\\])", "\\\\\\1", value)
}

generated_matches <- function(artifact) {
  basename_pattern <- paste0("^", escape_regex(basename(artifact)), "$")
  roots <- c(
    file.path("build", specification),
    file.path("build", "work", specification)
  )
  hits <- unlist(lapply(roots[dir.exists(roots)], function(root) {
    list.files(
      root, recursive = TRUE, full.names = TRUE,
      pattern = basename_pattern
    )
  }), use.names = FALSE)
  unique(hits)
}

artifact_lower <- tolower(registry$artifact)

owner <- rep("manual_or_static", nrow(registry))
owner[grepl("balance|attrition", artifact_lower)] <- "balance"
owner[grepl(
  "(^|/)rf-tables/|het-tes|takeup-np-rf|sms-te|dydc-derivative|dist-ri",
  artifact_lower
)] <- "reduced_form"
owner[grepl(
  "optim|plot-scaled|comp-dist-plot|ramsey-control|counterfactual-assignment",
  artifact_lower
)] <- "optimal_policy"
owner[grepl(
  "struct|wtp-summ|fob-beliefs|private-signal|fob-disagg|density_table",
  artifact_lower
) | grepl(
  "delta-with|sm-decomp|reported-social-perception|(^|/)dist-plot",
  artifact_lower
) | grepl(
  "w-distributions-robustness|gauss-mix|b-changes",
  artifact_lower
)] <- "structural_postprocess"
owner[grepl(
  "optim|plot-scaled|comp-dist-plot|ramsey-control",
  artifact_lower
)] <- "optimal_policy"
owner[grepl(
  "counterfactual-assignment|all-clusters-map",
  artifact_lower
)] <- "design_simulation"
owner[grepl(
  "signal_bracelet|calendar\\.png|sensitization|flyer_",
  artifact_lower
)] <- "versioned_source_asset"

generated_paths <- lapply(registry$artifact, generated_matches)
generated_by_build <- lengths(generated_paths) > 0L

coverage <- rep("manual_or_existing_only", nrow(registry))
coverage[owner == "versioned_source_asset" & registry$status == "resolved"] <-
  "covered_versioned_source"
coverage[owner %in% c("reduced_form", "balance") & generated_by_build] <-
  "covered_make_target"
coverage[owner %in% c("reduced_form", "balance") & !generated_by_build] <-
  "gap_expected_make_output"
coverage[owner == "structural_postprocess"] <-
  "gap_after_compact_gq"
coverage[owner == "optimal_policy"] <-
  "external_policy_workflow_not_make"
coverage[owner == "design_simulation"] <-
  "external_design_workflow_not_make"
coverage[registry$status == "missing"] <- "missing_source"

note <- rep("Staged from an existing source; no generator is tracked.", nrow(registry))
note[coverage == "covered_versioned_source"] <-
  "Versioned/static source asset; generation is not required."
note[coverage == "covered_make_target"] <-
  "Regenerated in the isolated Make/targets build."
note[coverage == "gap_expected_make_output"] <-
  "Owned by an existing Make analysis but not found among its isolated outputs."
note[coverage == "gap_after_compact_gq"] <-
  paste(
    "The pipeline creates compact-GQ intermediates but does not yet render",
    "this paper table or figure."
  )
note[coverage == "external_policy_workflow_not_make"] <-
  paste(
    "Requires the Gurobi/cluster policy workflow; make paper currently copies",
    "an existing artifact."
  )
note[coverage == "external_design_workflow_not_make"] <-
  paste(
    "Requires the study-design simulation/map workflow; make paper currently",
    "copies an existing artifact."
  )
note[coverage == "missing_source"] <-
  "The artifact is neither generated nor available in this repository snapshot."

registry$workflow_owner <- owner
registry$pipeline_coverage <- coverage
registry$generated_path <- vapply(
  generated_paths,
  function(paths) if (length(paths)) paths[[1L]] else NA_character_,
  character(1)
)
registry$coverage_note <- note

output <- "build/manifest/paper-pipeline-coverage.csv"
write.csv(registry, output, row.names = FALSE)

summary <- as.data.frame(table(
  workflow_owner = registry$workflow_owner,
  pipeline_coverage = registry$pipeline_coverage
))
summary <- summary[summary$Freq > 0L, ]
write.csv(
  summary,
  "build/manifest/paper-pipeline-coverage-summary.csv",
  row.names = FALSE
)

cat("Wrote", nrow(registry), "active dependencies to", output, "\n")
print(summary, row.names = FALSE)
