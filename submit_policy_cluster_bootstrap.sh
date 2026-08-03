#!/usr/bin/env bash

set -euo pipefail

stage_script=${STAGE_SCRIPT:-slurm_policy_cluster_bootstrap.sh}
common_export=${COMMON_EXPORT:-ALL}
num_replicates=${NUM_REPLICATES:-210}

prepare_job=$(sbatch --parsable \
  --export="${common_export},STAGE=prepare,NUM_REPLICATES=${num_replicates}" "${stage_script}")
predict_job=$(sbatch --parsable --dependency="afterok:${prepare_job}" \
  --export="${common_export},STAGE=predict,NUM_REPLICATES=${num_replicates}" "${stage_script}")
optimize_job=$(sbatch --parsable --dependency="afterok:${predict_job}" --array=1-5 \
  --export="${common_export},STAGE=optimize,NUM_REPLICATES=${num_replicates}" "${stage_script}")
summary_job=$(sbatch --parsable --dependency="afterok:${optimize_job}" \
  --export="${common_export},STAGE=summarize,NUM_REPLICATES=${num_replicates}" "${stage_script}")

echo "prepare=${prepare_job} predict=${predict_job} optimize=${optimize_job} summarize=${summary_job}"
