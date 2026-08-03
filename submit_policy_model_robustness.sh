#!/usr/bin/env bash

set -euo pipefail

models=(
  correct-observability second-order-observability grouped-lambda arm-lambda
  student-t5 cluster-shock
)
script=${STAGE_SCRIPT:-slurm_policy_model_robustness.sh}

for model in "${models[@]}"; do
  prepare=$(sbatch --parsable --export="ALL,STAGE=prepare,MODEL_ID=${model}" "${script}")
  predict=$(sbatch --parsable --dependency="afterok:${prepare}" \
    --export="ALL,STAGE=predict,MODEL_ID=${model}" "${script}")
  optimize=$(sbatch --parsable --mem=4G --dependency="afterok:${predict}" --array=1-5 \
    --export="ALL,STAGE=optimize,MODEL_ID=${model}" "${script}")
  summarize=$(sbatch --parsable --mem=4G --dependency="afterok:${optimize}" \
    --export="ALL,STAGE=summarize,MODEL_ID=${model}" "${script}")
  echo "${model}: prepare=${prepare} predict=${predict} optimize=${optimize} summarize=${summarize}"
done
