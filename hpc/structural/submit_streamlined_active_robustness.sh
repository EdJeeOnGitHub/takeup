#!/usr/bin/env bash
set -euo pipefail

fit_script=${FIT_SCRIPT:-hpc/structural/slurm_streamlined_active_robustness.sh}
gq_script=${GQ_SCRIPT:-hpc/structural/slurm_streamlined_active_robustness_gq.sh}
max_fit_jobs=${MAX_FIT_JOBS:-8}

fit_job=$(sbatch --parsable --array=1-16%${max_fit_jobs} "${fit_script}")
gq_job=$(sbatch --parsable --dependency=afterok:${fit_job} --array=1-4 "${gq_script}")

printf '%-18s %s\n' fit "${fit_job}" compact_gq_1250 "${gq_job}"
