#!/usr/bin/env bash
set -euo pipefail

worker=${WORKER:-slurm_main_core_cluster_bootstrap_999.sh}
batch_size=${BATCH_SIZE:-333}
max_concurrent=${MAX_CONCURRENT:-32}
weights_job=$(sbatch --parsable --export=ALL,STAGE=weights "${worker}")
baseline_job=$(sbatch --parsable --dependency="afterok:${weights_job}" \
  --export=ALL,STAGE=baseline "${worker}")
first_end=$((batch_size < 999 ? batch_size : 999))
bootstrap_job=$(sbatch --parsable --dependency="afterok:${baseline_job}" \
  --array="1-${first_end}%${max_concurrent}" \
  --export=ALL,STAGE=bootstrap "${worker}")
launcher_job=$(sbatch --parsable --dependency="afterok:${bootstrap_job}" \
  --export=ALL,STAGE=launch,NEXT_START=$((first_end + 1)),BATCH_SIZE="${batch_size}",MAX_CONCURRENT="${max_concurrent}" \
  "${worker}")
printf '%-12s %s\n' weights "${weights_job}" baseline "${baseline_job}" \
  bootstrap_1 "${bootstrap_job}" launcher "${launcher_job}"
