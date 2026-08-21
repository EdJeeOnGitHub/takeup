#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (!length(hit)) default else substring(hit[[1L]], nchar(name) + 2L)
}
specification <- value("--spec", Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned"))
source("R/distance/spec.R")
specification <- takeup_distance_spec(specification)
output_root <- value("--output-root", file.path("build", "candidate-hpc", specification, "export"))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]]
tracked_status <- system2(
  "git", c("status", "--porcelain", "--untracked-files=no"), stdout = TRUE
)
allow_dirty <- "--allow-dirty" %in% args
if (length(tracked_status) && !allow_dirty) {
  stop(
    "Refusing to export external jobs from uncommitted tracked code. Commit ",
    "the implementation first, or use --allow-dirty only for a smoke test.",
    call. = FALSE
  )
}
sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  sub("[[:space:]].*$", "", system2("sha256sum", path, stdout = TRUE)[[1L]])
}
inputs <- c(
  "data/rct_targetable_schools_2.0.rds",
  "data/takeup_processed_cluster_strat.rds",
  "data/stan_analysis_data/dist_fit104.RData",
  "stan_models/takeup_struct_main_core.stan",
  "stan_models/takeup_struct_main_core_compact_gq.stan"
)
input_manifest <- data.frame(path = inputs, exists = file.exists(inputs),
                             sha256 = vapply(inputs, sha256, character(1)))
utils::write.csv(input_manifest, file.path(output_root, "input-manifest.csv"),
                 row.names = FALSE, na = "")

workflows <- data.frame(
  workflow_id = c(
    "baseline", "cluster-shock", "cluster-weight-999", "prior-grid",
    "student-t", "finite-mixture", "lambda", "observability-reporting",
    "individual-distance", "no-outliers", "policy-parameters",
    "policy-model-inputs"
  ),
  required = TRUE,
  command = c(
    "make structural-fit DISTANCE_SPEC=assigned",
    "DISTANCE_DEFINITION=assigned USE_CORE_CLUSTER_SHOCK=1 sbatch --array=1-4 hpc/structural/slurm_main_core.sh",
    "DISTANCE_DEFINITION=assigned bash hpc/structural/submit_main_core_cluster_bootstrap_999.sh",
    "DISTANCE_DEFINITION=assigned bash hpc/structural/submit_main_core_prior_grid.sh",
    "DISTANCE_DEFINITION=assigned sbatch hpc/structural/slurm_main_core_student_t.sh",
    "DISTANCE_DEFINITION=assigned bash hpc/structural/submit_main_core_finite_mixture_midway3.sh",
    "DISTANCE_DEFINITION=assigned sbatch hpc/structural/slurm_main_core_lambda.sh",
    "DISTANCE_DEFINITION=assigned sbatch hpc/structural/slurm_main_core_observability_gq.sh",
    "DISTANCE_DEFINITION=assigned MODEL=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP sbatch --array=1-4 hpc/structural/slurm_main_core.sh",
    "DISTANCE_DEFINITION=assigned MODEL=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS sbatch --array=1-4 hpc/structural/slurm_main_core.sh",
    "Rscript scripts/policy/prepare-cluster-bootstrap.R --distance-definition=assigned --num-replicates=999",
    "Rscript scripts/policy/prepare-model-robustness.R and predict-model-robustness.R for every maintained alternative"
  ),
  stringsAsFactors = FALSE
)
expected <- c(
  baseline = "tables/struct-overall-te-table.tex;tables/private-signal-te-table.tex;tables/fob-beliefs-table.tex;tables/wtp-summ-table.tex",
  `cluster-shock` = "appendix/structural-robustness/tables/main-core-cluster-robustness.tex",
  `cluster-weight-999` = paste(c(
    "appendix/structural-robustness/tables/main-core-exponential-cluster-weight-overall-te-table.tex",
    "appendix/structural-robustness/tables/main-core-exponential-cluster-weight-multiplier-contrasts.tex",
    "appendix/structural-robustness/tables/main-core-exponential-cluster-weight-finite-segments.tex"
  ), collapse = ";"),
  `prior-grid` = "appendix/structural-robustness/tables/main-core-prior-grid.tex",
  `student-t` = "appendix/structural-robustness/tables/main-core-student-t-multipliers.tex",
  `finite-mixture` = "appendix/structural-robustness/tables/main-core-finite-mixture-multipliers.tex;appendix/structural-robustness/figures/finite-mixture-density.pdf",
  lambda = "appendix/structural-robustness/tables/main-core-lambda-nested-ates.tex;appendix/structural-robustness/tables/main-core-lambda-nested-multipliers.tex",
  `observability-reporting` = "appendix/structural-robustness/tables/main-core-observability-multipliers.tex;appendix/structural-robustness/tables/main-core-tight-asymmetric-report-multipliers.tex;appendix/structural-robustness/tables/main-core-report-distance-contrasts.tex",
  `individual-distance` = "tables/indiv-dist-community-fp-indiv-vis-robust-struct-overall-te-table.tex;tables/indiv-dist-indiv-fp-robust-struct-overall-te-table.tex",
  `no-outliers` = "tables/struct-robustness-nooutliers-overall-te-table.tex",
  `policy-parameters` = "policy/policy-bootstrap-parameters.csv;policy/policy-bootstrap-manifest.csv",
  `policy-model-inputs` = "policy/model-robustness/correct-observability/policy-model-summary.csv;policy/model-robustness/second-order-observability/policy-model-summary.csv;policy/model-robustness/grouped-lambda/policy-model-summary.csv;policy/model-robustness/arm-lambda/policy-model-summary.csv;policy/model-robustness/student-t5/policy-model-summary.csv;policy/model-robustness/cluster-shock/policy-model-summary.csv"
)
workflows$expected_artifacts <- unname(expected[workflows$workflow_id])
utils::write.csv(workflows, file.path(output_root, "job-manifest.csv"), row.names = FALSE)
manifest <- data.frame(
  schema_version = 1L, distance_definition = specification,
  git_commit = commit, generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  workflow_count = nrow(workflows), tracked_worktree_dirty = length(tracked_status) > 0L
)
utils::write.csv(manifest, file.path(output_root, "candidate-hpc-manifest.csv"), row.names = FALSE)
writeLines(c(
  "# Candidate HPC bundle", "",
  "Run every command in `job-manifest.csv` from the recorded Git commit.",
  "The completed bundle must contain an unchanged candidate-hpc-manifest.csv,",
  "an artifact-manifest.csv with columns path,sha256,workflow_id,distance_definition,",
  "and an artifacts/ directory containing those relative paths."
), file.path(output_root, "README.md"))
cat(output_root, "\n")
