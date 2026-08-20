#!/usr/bin/env bash
set -euo pipefail

script=${SCRIPT:-slurm_main_core_prior_grid.sh}
max_sample_jobs=${MAX_SAMPLE_JOBS:-16}
max_gq_jobs=${MAX_GQ_JOBS:-8}
resume_from_compiled=${RESUME_FROM_COMPILED:-0}

if [[ "${resume_from_compiled}" == 1 ]]; then
  prepare_job=skipped
  compile_job=skipped
  sample_job=$(sbatch --parsable --array=1-52%${max_sample_jobs} \
    --export=ALL,STAGE=sample "${script}")
else
  prepare_job=$(sbatch --parsable --export=ALL,STAGE=prepare "${script}")
  compile_job=$(sbatch --parsable --dependency=afterok:${prepare_job} \
    --export=ALL,STAGE=compile "${script}")
  sample_job=$(sbatch --parsable --dependency=afterok:${compile_job} \
    --array=1-52%${max_sample_jobs} --export=ALL,STAGE=sample "${script}")
fi
gq_job=$(sbatch --parsable --dependency=afterok:${sample_job} \
  --array=1-13%${max_gq_jobs} --export=ALL,STAGE=gq "${script}")
summary_job=$(sbatch --parsable --dependency=afterok:${gq_job} \
  --export=ALL,STAGE=summarize "${script}")

printf '%-18s %s\n' \
  prepare "${prepare_job}" \
  compile "${compile_job}" \
  sample "${sample_job}" \
  compact_gq "${gq_job}" \
  summarize "${summary_job}"
