#!/usr/bin/env bash
set -euo pipefail

project_root=${PROJECT_ROOT:-$(pwd)}
worker=${WORKER:-${project_root}/slurm_main_core_cluster_bootstrap_999.sh}
batch_size=${BATCH_SIZE:-333}
max_concurrent=${MAX_CONCURRENT:-32}
total_replicates=${NUM_REPLICATES:-999}

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
export WORKER_SCRIPT="${worker}"
export TOTAL_REPLICATES="${total_replicates}"
export SLURM_PARTITION_NAME="${slurm_partition}"
export SLURM_ACCOUNT_NAME="${slurm_account}"

weights_job=$(sbatch --parsable "${sbatch_site[@]}" \
  --export=ALL,STAGE=weights "${worker}")
baseline_job=$(sbatch --parsable "${sbatch_site[@]}" \
  --dependency="afterok:${weights_job}" \
  --export=ALL,STAGE=baseline "${worker}")
first_end=$((batch_size < total_replicates ? batch_size : total_replicates))
bootstrap_job=$(sbatch --parsable "${sbatch_site[@]}" \
  --dependency="afterok:${baseline_job}" \
  --array="1-${first_end}%${max_concurrent}" \
  --export=ALL,STAGE=bootstrap "${worker}")
launcher_job=$(sbatch --parsable "${sbatch_site[@]}" \
  --dependency="afterok:${bootstrap_job}" \
  --export=ALL,STAGE=launch,NEXT_START=$((first_end + 1)),BATCH_SIZE="${batch_size}",MAX_CONCURRENT="${max_concurrent}" \
  "${worker}")
printf '%-12s %s\n' weights "${weights_job}" baseline "${baseline_job}" \
  bootstrap_1 "${bootstrap_job}" launcher "${launcher_job}"
