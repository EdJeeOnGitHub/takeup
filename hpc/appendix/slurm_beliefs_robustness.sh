#!/usr/bin/env bash

#SBATCH --partition=broadwl
#SBATCH --job-name=takeup-beliefs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-10:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-beliefs-robustness-%A_%a.log
#SBATCH --error=temp/log/takeup-beliefs-robustness-%A_%a.log
#SBATCH --export=IN_SLURM=1
#SBATCH --array=0-2

set -euo pipefail

LATEST_VERSION=106
VERSION=${1:-$LATEST_VERSION}
CMDSTAN_ARGS="--cmdstanr"
SLURM_INOUT_DIR="/project/akaring/takeup-data/data/stan_analysis_data"
ITER=${ITER:-400}

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."
  module load -f midway2 gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0

  OUTPUT_ARGS="--output-path=${SLURM_INOUT_DIR}"
  POSTPROCESS_INOUT_ARGS="--input-path=${SLURM_INOUT_DIR} --output-path=${SLURM_INOUT_DIR}"
  CORES=$SLURM_CPUS_PER_TASK
else
  OUTPUT_ARGS="--output-path=data/stan_analysis_data"
  POSTPROCESS_INOUT_ARGS=
  CORES=8
fi

models=(
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB_MISSING_MARGINALIZED"
)

model=${models[$SLURM_ARRAY_TASK_ID]}
STAN_THREADS=$((${CORES} / 4))

source scripts/structural/run-postprocess.sh

Rscript --no-save \
  --no-restore \
  --verbose \
  scripts/structural/run-model.R takeup fit \
  --models=${model} \
  ${CMDSTAN_ARGS} \
  ${OUTPUT_ARGS} \
  --update-output \
  --threads=${STAN_THREADS} \
  --outputname=dist_fit${VERSION} \
  --num-mix-groups=1 \
  --chains=4 \
  --iter=${ITER} \
  --sequential > temp/log/output-${model}-fit${VERSION}.txt 2>&1

postprocess_struct_models "${model}" "${VERSION}" "${POSTPROCESS_INOUT_ARGS}"
