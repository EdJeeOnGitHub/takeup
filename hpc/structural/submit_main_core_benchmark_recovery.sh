#!/usr/bin/env bash
set -euo pipefail

MAX_JOBS=${MAX_JOBS:-20}
: "${BASELINE_CSV_DIR:?Set BASELINE_CSV_DIR to the latest assigned-distance benchmark chain directory}"
generate_job=$(sbatch --parsable --export=ALL,STAGE=generate \
  hpc/structural/slurm_main_core_benchmark_recovery.sh)
compile_job=$(sbatch --parsable --dependency=afterok:${generate_job} \
  --export=ALL,STAGE=compile hpc/structural/slurm_main_core_benchmark_recovery.sh)
fit_job=$(sbatch --parsable --dependency=afterok:${compile_job} \
  --array=1-50%${MAX_JOBS} --export=ALL,STAGE=fit \
  hpc/structural/slurm_main_core_benchmark_recovery.sh)
summary_job=$(sbatch --parsable --dependency=afterok:${fit_job} \
  --export=ALL,STAGE=summarize hpc/structural/slurm_main_core_benchmark_recovery.sh)

printf '%-18s %s\n' generate "${generate_job}" compile "${compile_job}" \
  fits "${fit_job}" summarize "${summary_job}"
