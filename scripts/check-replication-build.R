#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- sub("^--spec=", "", args[grepl("^--spec=", args)])
if (!length(value)) value <- Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned")
source("R/distance-spec.R")
specification <- takeup_distance_spec(value[[1L]])
crosswalk <- takeup_distance_crosswalk()
analysis_crosswalk <- dplyr::filter(crosswalk, .data$in_main_analysis)
stopifnot(nrow(analysis_crosswalk) == 144L,
          sum(analysis_crosswalk$switched) == 26L)

registry_path <- "build/manifest/paper-artifacts.csv"
if (!file.exists(registry_path)) stop("Run the artifact registry first.")
registry <- read.csv(registry_path, stringsAsFactors = FALSE)
coverage_path <- "build/manifest/paper-pipeline-coverage.csv"
if (!file.exists(coverage_path)) stop("Run the paper pipeline audit first.")
coverage <- read.csv(coverage_path, stringsAsFactors = FALSE)

required_scripts <- c(
  "scratch/sample-slim-individual-fp.R",
  "stan_models/takeup_struct_main_core.stan",
  "scratch/generate-main-core-compact-gq.R",
  "stan_models/takeup_struct_main_core_compact_gq.stan"
)
if (!all(file.exists(required_scripts))) stop("Latest structural workflow is incomplete.")

dir.create(file.path("build", specification, "audit"), recursive = TRUE,
           showWarnings = FALSE)
checks <- data.frame(
  check = c("clusters", "switched_clusters", "paper_dependencies",
            "unresolved_dependencies", "latest_structural_files",
            "default_covered_dependencies", "default_uncovered_dependencies",
            "structural_full_render_gaps", "external_policy_dependencies"),
  value = c(nrow(analysis_crosswalk), sum(analysis_crosswalk$switched), nrow(registry),
            sum(registry$status == "missing"), sum(file.exists(required_scripts)),
            sum(grepl("^covered_", coverage$pipeline_coverage)),
            sum(!grepl("^covered_", coverage$pipeline_coverage)),
            sum(coverage$full_regeneration_coverage == "render_target_gap"),
            sum(coverage$full_regeneration_coverage ==
                  "external_generator_available"))
)
path <- file.path("build", specification, "audit", "replication-checks.csv")
write.csv(checks, path, row.names = FALSE)
print(checks, row.names = FALSE)
if (sum(!grepl("^covered_", coverage$pipeline_coverage)) > 0L) {
  stop("The default paper contract has uncovered dependencies; see ", coverage_path)
}
