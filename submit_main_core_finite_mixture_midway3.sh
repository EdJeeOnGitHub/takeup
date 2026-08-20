#!/usr/bin/env bash
set -euo pipefail

script=${SCRIPT:-slurm_main_core_finite_mixture_midway3.sh}
compile_job=$(sbatch --parsable --export=ALL,STAGE=compile "${script}")
sample_job=$(sbatch --parsable --dependency=afterok:${compile_job} --array=1-4 \
  --export=ALL,STAGE=sample "${script}")
gq_job=$(sbatch --parsable --dependency=afterok:${sample_job} \
  --export=ALL,STAGE=gq "${script}")
summary_job=$(sbatch --parsable --dependency=afterok:${gq_job} \
  --export=ALL,STAGE=summarize "${script}")

printf '%-20s %s\n' \
  compile "${compile_job}" \
  sample "${sample_job}" \
  gq "${gq_job}" \
  summary "${summary_job}"
