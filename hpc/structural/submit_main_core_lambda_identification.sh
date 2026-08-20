#!/usr/bin/env bash

set -euo pipefail

PHASE=${1:-core}
MAX_PROFILE_JOBS=${MAX_PROFILE_JOBS:-20}
MAX_RECOVERY_JOBS=${MAX_RECOVERY_JOBS:-40}
MAX_HMC_JOBS=${MAX_HMC_JOBS:-20}

case "${PHASE}" in
  core)
    prepare_job=$(sbatch --parsable --export=ALL,STAGE=prepare \
      hpc/structural/slurm_main_core_lambda.sh)
    reuse_job=$(sbatch --parsable --dependency=afterok:${prepare_job} \
      --export=ALL,STAGE=reuse_common hpc/structural/slurm_main_core_lambda.sh)
    # Compile both Stan executables once before any array can touch them.
    gate_job=$(sbatch --parsable --dependency=afterok:${reuse_job} --array=1 \
      --export=ALL,STAGE=mode hpc/structural/slurm_main_core_lambda.sh)
    mode_job=$(sbatch --parsable --dependency=afterok:${gate_job} --array=2-9 \
      --export=ALL,STAGE=mode hpc/structural/slurm_main_core_lambda.sh)
    sample_job=$(sbatch --parsable \
      --dependency=afterok:${gate_job}:${mode_job} --array=5-36 \
      --export=ALL,STAGE=sample hpc/structural/slurm_main_core_lambda.sh)
    gq_job=$(sbatch --parsable \
      --dependency=afterany:${sample_job}:${reuse_job} --array=1-9 \
      --export=ALL,STAGE=gq hpc/structural/slurm_main_core_lambda.sh)
    grid_summary_job=$(sbatch --parsable --dependency=afterany:${gq_job} \
      --export=ALL,STAGE=summarize hpc/structural/slurm_main_core_lambda.sh)
    profile_prepare_job=$(sbatch --parsable \
      --dependency=afterok:${gate_job}:${mode_job} \
      --export=ALL,STAGE=prepare hpc/structural/slurm_main_core_lambda_profile.sh)
    profile_job=$(sbatch --parsable \
      --dependency=afterok:${profile_prepare_job} \
      --array=1-61%${MAX_PROFILE_JOBS} --export=ALL,STAGE=profile \
      hpc/structural/slurm_main_core_lambda_profile.sh)
    profile_summary_job=$(sbatch --parsable --dependency=afterany:${profile_job} \
      --export=ALL,STAGE=summarize hpc/structural/slurm_main_core_lambda_profile.sh)
    printf '%-24s %s\n' \
      prepare "${prepare_job}" reuse_common "${reuse_job}" \
      compile_mode_gate "${gate_job}" modes "${mode_job}" \
      production_samples "${sample_job}" compact_gq "${gq_job}" \
      grid_summary "${grid_summary_job}" \
      profile_prepare "${profile_prepare_job}" profile "${profile_job}" \
      profile_summary "${profile_summary_job}"
    ;;
  recovery-generate)
    sbatch --export=ALL,STAGE=generate hpc/structural/slurm_main_core_lambda_recovery.sh
    ;;
  recovery-fit-a)
    sbatch --array=1-375%${MAX_RECOVERY_JOBS} --export=ALL,STAGE=fit \
      hpc/structural/slurm_main_core_lambda_recovery.sh
    ;;
  recovery-fit-b)
    sbatch --array=376-750%${MAX_RECOVERY_JOBS} --export=ALL,STAGE=fit \
      hpc/structural/slurm_main_core_lambda_recovery.sh
    ;;
  recovery-hmc)
    sbatch --array=1-120%${MAX_HMC_JOBS} --export=ALL,STAGE=hmc \
      hpc/structural/slurm_main_core_lambda_recovery.sh
    ;;
  recovery-finish)
    hmc_gq_job=$(sbatch --parsable --array=1-30 \
      --export=ALL,STAGE=hmc_gq hpc/structural/slurm_main_core_lambda_recovery.sh)
    summary_job=$(sbatch --parsable --dependency=afterany:${hmc_gq_job} \
      --export=ALL,STAGE=summarize hpc/structural/slurm_main_core_lambda_recovery.sh)
    assemble_job=$(sbatch --parsable --dependency=afterok:${summary_job} \
      --export=ALL,STAGE=assemble hpc/structural/slurm_main_core_lambda_recovery.sh)
    printf '%-24s %s\n' hmc_gq "${hmc_gq_job}" \
      recovery_summary "${summary_job}" paper_package "${assemble_job}"
    ;;
  *)
    echo "Unknown phase: ${PHASE}" >&2
    echo "Use core, recovery-generate, recovery-fit-a, recovery-fit-b, recovery-hmc, or recovery-finish." >&2
    exit 2
    ;;
esac
