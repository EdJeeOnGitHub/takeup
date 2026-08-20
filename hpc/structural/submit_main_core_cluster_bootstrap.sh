#!/usr/bin/env bash

set -euo pipefail

stage_script=${STAGE_SCRIPT:-hpc/structural/slurm_main_core_cluster_bootstrap.sh}
common_export=${COMMON_EXPORT:-ALL}

prepare_job=$(sbatch --parsable \
  --export="${common_export},STAGE=prepare" "${stage_script}")
baseline_job=$(sbatch --parsable --dependency="afterok:${prepare_job}" \
  --export="${common_export},STAGE=baseline" "${stage_script}")
bootstrap_job=$(sbatch --parsable --dependency="afterok:${baseline_job}" \
  --array=1-200%32 --export="${common_export},STAGE=bootstrap" \
  "${stage_script}")
summary_job=$(sbatch --parsable --dependency="afterok:${bootstrap_job}" \
  --export="${common_export},STAGE=summarize" "${stage_script}")

echo "prepare=${prepare_job} baseline=${baseline_job} bootstrap=${bootstrap_job} summarize=${summary_job}"
