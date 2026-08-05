#!/usr/bin/env bash
set -euo pipefail

INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-asym-input}
OUTPUT_ROOT=${OUTPUT_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors}
OUTPUT_PATH=${OUTPUT_ROOT}/design-pooled
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS}

sample_job=$(sbatch --parsable --array=1-4 \
  --export="ALL,MODEL=${MODEL},INPUT_PATH=${INPUT_PATH},OUTPUT_PATH=${OUTPUT_PATH},STAN_PATH=stan_models,CORE_OBSERVATION_MODEL=1,CORE_RECOGNITION_STRUCTURE=0,CORE_REPORT_STRUCTURE=0,CORE_REPORT_ARM_DIST_HIERARCHICAL=2,CORE_REPORT_ARM_DIST_PRIOR_SCALE=0.25,ITER_WARMUP=800,ITER_SAMPLING=800,THREADS_PER_CHAIN=8,ADAPT_DELTA=0.99,MAX_TREEDEPTH=12,METRIC=diag_e,REFRESH=25" \
  slurm_main_core.sh)
gq_job=$(sbatch --parsable --dependency="afterok:${sample_job}" --array=1-4 \
  --export="ALL,SPECIFICATION=design-pooled,FIT_ROOT=${OUTPUT_ROOT},INPUT_PATH=${INPUT_PATH}" \
  slurm_main_core_observability_gq.sh)
summary_job=$(sbatch --parsable --dependency="afterok:${gq_job}" \
  --export="ALL,ROOT=${OUTPUT_ROOT}" slurm_main_core_report_distance_summary.sh)
printf '%-12s %s\n' sample "${sample_job}" gq "${gq_job}" summary "${summary_job}"
