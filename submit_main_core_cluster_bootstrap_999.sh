#!/usr/bin/env bash
set -euo pipefail

worker=${WORKER:-slurm_main_core_cluster_bootstrap_999.sh}
weights_job=$(sbatch --parsable --export=ALL,STAGE=weights "${worker}")
baseline_job=$(sbatch --parsable --dependency="afterok:${weights_job}" \
  --export=ALL,STAGE=baseline "${worker}")
bootstrap_job=$(sbatch --parsable --dependency="afterok:${baseline_job}" \
  --array=1-999%32 --export=ALL,STAGE=bootstrap "${worker}")
summary_job=$(sbatch --parsable --dependency="afterok:${bootstrap_job}" \
  --export=ALL,STAGE=summarize "${worker}")
printf '%-12s %s\n' weights "${weights_job}" baseline "${baseline_job}" \
  bootstrap "${bootstrap_job}" summary "${summary_job}"
