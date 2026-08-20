#!/usr/bin/env bash
set -euo pipefail

script=${SCRIPT:-slurm_main_core_prior_grid.sh}
output_path=${OUTPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-prior-grid}
audit=${RERUN_AUDIT:-${output_path}/prior-grid-needs-rerun.csv}
max_sample_jobs=${MAX_SAMPLE_JOBS:-16}
max_gq_jobs=${MAX_GQ_JOBS:-8}

module load -f R/4.2.0

if [[ ! -f "${audit}" ]]; then
  echo "Missing initial-run audit: ${audit}" >&2
  exit 2
fi
n_specs=$(Rscript --no-save --no-restore --no-init-file -e \
  "x <- read.csv('${audit}'); cat(nrow(x))")
if [[ "${n_specs}" -eq 0 ]]; then
  echo "No specifications are flagged for an 800/800 rerun."
  exit 0
fi

sample_job=$(sbatch --parsable --array=1-$((4 * n_specs))%${max_sample_jobs} \
  --export=ALL,STAGE=sample-rerun,ITER_WARMUP=800,ITER_SAMPLING=800 \
  "${script}")
gq_job=$(sbatch --parsable --dependency=afterok:${sample_job} \
  --array=1-${n_specs}%${max_gq_jobs} \
  --export=ALL,STAGE=gq-rerun "${script}")
summary_job=$(sbatch --parsable --dependency=afterok:${gq_job} \
  --export=ALL,STAGE=summarize "${script}")

printf '%-18s %s\n' \
  rerun_sample "${sample_job}" \
  rerun_compact_gq "${gq_job}" \
  final_summary "${summary_job}"
