#!/usr/bin/env bash
set -euo pipefail
mode_job=$(sbatch --parsable --export=ALL,STAGE=mode slurm_main_core_student_t.sh)
sample_job=$(sbatch --parsable --dependency=afterok:${mode_job} --array=1-4 \
  --export=ALL,STAGE=sample slurm_main_core_student_t.sh)
gq_job=$(sbatch --parsable --dependency=afterok:${sample_job} \
  --export=ALL,STAGE=gq slurm_main_core_student_t.sh)
summary_job=$(sbatch --parsable --dependency=afterok:${gq_job} \
  --export=ALL,STAGE=summarize slurm_main_core_student_t.sh)
printf '%-20s %s\n' mode "${mode_job}" sample "${sample_job}" \
  gq "${gq_job}" summary "${summary_job}"
