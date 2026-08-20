#!/usr/bin/env bash

set -euo pipefail

project_root=${PROJECT_ROOT:-$(pwd)}
stage_script=${STAGE_SCRIPT:-${project_root}/hpc/policy/slurm_policy_cluster_bootstrap.sh}
common_export=${COMMON_EXPORT:-ALL}
num_replicates=${NUM_REPLICATES:-999}

cluster_name=${MIDWAY_CLUSTER:-$(
  scontrol show config 2>/dev/null | awk '/^ClusterName/ {print $3}'
)}
if [[ "${cluster_name}" == "midway3" ]]; then
  slurm_partition=${SLURM_PARTITION_NAME:-caslake}
  slurm_account=${SLURM_ACCOUNT_NAME:-pi-akaring}
else
  slurm_partition=${SLURM_PARTITION_NAME:-broadwl}
  slurm_account=${SLURM_ACCOUNT_NAME:-}
fi
sbatch_site=(--partition="${slurm_partition}")
if [[ -n "${slurm_account}" ]]; then
  sbatch_site+=(--account="${slurm_account}")
fi

export PROJECT_ROOT="${project_root}"
export SLURM_PARTITION_NAME="${slurm_partition}"
export SLURM_ACCOUNT_NAME="${slurm_account}"

prepare_job=$(sbatch --parsable "${sbatch_site[@]}" \
  --export="${common_export},STAGE=prepare,NUM_REPLICATES=${num_replicates}" "${stage_script}")
predict_job=$(sbatch --parsable "${sbatch_site[@]}" --dependency="afterok:${prepare_job}" \
  --export="${common_export},STAGE=predict,NUM_REPLICATES=${num_replicates}" "${stage_script}")
optimize_job=$(sbatch --parsable "${sbatch_site[@]}" --dependency="afterok:${predict_job}" --array=1-5 \
  --export="${common_export},STAGE=optimize,NUM_REPLICATES=${num_replicates}" "${stage_script}")
summary_job=$(sbatch --parsable "${sbatch_site[@]}" --dependency="afterok:${optimize_job}" \
  --export="${common_export},STAGE=summarize,NUM_REPLICATES=${num_replicates}" "${stage_script}")
population_job=$(sbatch --parsable "${sbatch_site[@]}" --dependency="afterok:${summary_job}" \
  --export="${common_export},STAGE=population,NUM_REPLICATES=${num_replicates}" "${stage_script}")

echo "prepare=${prepare_job} predict=${predict_job} optimize=${optimize_job} summarize=${summary_job} population=${population_job}"
