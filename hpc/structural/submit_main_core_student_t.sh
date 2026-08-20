#!/usr/bin/env bash
set -euo pipefail
mode_job=$(sbatch --parsable --export=ALL,STAGE=mode hpc/structural/slurm_main_core_student_t.sh)
sample_job=$(sbatch --parsable --dependency=afterok:${mode_job} --array=1-4 \
  --export=ALL,STAGE=sample hpc/structural/slurm_main_core_student_t.sh)
gq_job=$(sbatch --parsable --dependency=afterok:${sample_job} \
  --export=ALL,STAGE=gq hpc/structural/slurm_main_core_student_t.sh)
summary_job=$(sbatch --parsable --dependency=afterok:${gq_job} \
  --export=ALL,STAGE=summarize hpc/structural/slurm_main_core_student_t.sh)
printf '%-20s %s\n' mode "${mode_job}" sample "${sample_job}" \
  gq "${gq_job}" summary "${summary_job}"
