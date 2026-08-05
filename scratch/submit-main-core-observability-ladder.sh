#!/usr/bin/env bash

# Submit the restricted asymmetric-observability ladder. F0 is the already
# completed unrestricted conditional fit and is intentionally not rerun here.
set -euo pipefail

INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-asym-input}
OUTPUT_ROOT=${OUTPUT_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS}
ITER_WARMUP=${ITER_WARMUP:-400}
ITER_SAMPLING=${ITER_SAMPLING:-400}
THREADS_PER_CHAIN=${THREADS_PER_CHAIN:-8}
ADAPT_DELTA=${ADAPT_DELTA:-0.99}
MAX_TREEDEPTH=${MAX_TREEDEPTH:-12}

# label:observation model:recognition structure:report structure
specifications=(
  f1:1:0:1
  f2:1:2:1
  f3:1:2:2
  u3:2:1:2
)

for specification in "${specifications[@]}"; do
  IFS=: read -r label observation recognition report <<< "${specification}"
  job_id=$(sbatch --parsable --array=1-4 \
    --export="ALL,MODEL=${MODEL},INPUT_PATH=${INPUT_PATH},OUTPUT_PATH=${OUTPUT_ROOT}/${label},STAN_PATH=stan_models,CORE_OBSERVATION_MODEL=${observation},CORE_RECOGNITION_STRUCTURE=${recognition},CORE_REPORT_STRUCTURE=${report},ITER_WARMUP=${ITER_WARMUP},ITER_SAMPLING=${ITER_SAMPLING},THREADS_PER_CHAIN=${THREADS_PER_CHAIN},ADAPT_DELTA=${ADAPT_DELTA},MAX_TREEDEPTH=${MAX_TREEDEPTH},METRIC=diag_e,REFRESH=25" \
    slurm_main_core.sh)
  printf '%s:%s\n' "${label}" "${job_id}"
done
