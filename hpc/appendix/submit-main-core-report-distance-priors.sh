#!/usr/bin/env bash

# Compare partial pooling against a tighter independent prior while preserving
# the original fully flexible conditional multinomial likelihood.
set -euo pipefail

INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-asym-input}
OUTPUT_ROOT=${OUTPUT_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS}
ITER_WARMUP=${ITER_WARMUP:-400}
ITER_SAMPLING=${ITER_SAMPLING:-400}

# label:hierarchical:prior scale
specifications=(hierarchical:1:0.25 tight:0:0.10)
for specification in "${specifications[@]}"; do
  IFS=: read -r label hierarchical prior_scale <<< "${specification}"
  job_id=$(sbatch --parsable --array=1-4 \
    --export="ALL,MODEL=${MODEL},INPUT_PATH=${INPUT_PATH},OUTPUT_PATH=${OUTPUT_ROOT}/${label},STAN_PATH=stan_models,CORE_OBSERVATION_MODEL=1,CORE_RECOGNITION_STRUCTURE=0,CORE_REPORT_STRUCTURE=0,CORE_REPORT_ARM_DIST_HIERARCHICAL=${hierarchical},CORE_REPORT_ARM_DIST_PRIOR_SCALE=${prior_scale},ITER_WARMUP=${ITER_WARMUP},ITER_SAMPLING=${ITER_SAMPLING},THREADS_PER_CHAIN=8,ADAPT_DELTA=0.99,MAX_TREEDEPTH=12,METRIC=diag_e,REFRESH=25" \
    hpc/structural/slurm_main_core.sh)
  printf '%s:%s\n' "${label}" "${job_id}"
done
