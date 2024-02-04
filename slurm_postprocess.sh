#!/bin/bash
#SBATCH --partition=bigmem2
#SBATCH --job-name=quick-takeup
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=0-10:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-%j.log
#SBATCH --error=temp/log/takeup-%j.log
#SBATCH --export=IN_SLURM=1

LATEST_VERSION=101
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided
SLURM_INOUT_DIR="/project/akaring/takeup-data/"

models=(
   "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
   "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_FOB"
   "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_FIXED_FOB"
#  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_SD_WTP_VAL"
#  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_MU_WTP_VAL"
#  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_WTP_SUBMODEL"
#  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_BELIEFS_SUBMODEL"
#  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_SUBMODELS"
)

echo "Version: $VERSION"

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."

  module load -f midway2 gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
  module load -f R/4.2.0
  IN_ARG="--input-path=${SLURM_INOUT_DIR}/data/stan_analysis_data"
  OUT_ARG="--output-path=${SLURM_INOUT_DIR}/temp-data"

  echo "Running with ${CORES} cores."
  echo "INOUT ARGS: ${IN_ARG} ${OUT_ARG}."
else
  IN_ARG="--input-path=data/stan_analysis_data"
  OUT_ARG="--output-path=temp-data"
  CORES=8
fi

for model in "${models[@]}"; do
  Rscript --no-save \
          --no-restore \
          --verbose \
          quick_ate_postprocess.R \
          ${VERSION} \
          --model=${model} \
	  ${IN_ARG} \
	  ${OUT_ARG} \
          1 2 3 4  &
done

wait
