#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --job-name=takeup-indiv
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-indiv-%j.log
#SBATCH --error=temp/log/takeup-indiv-%j.log
#SBATCH --export=IN_SLURM=1

VERSION=105
ITER=400
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP"
SLURM_INOUT_DIR="/project/akaring/takeup-data/data/stan_analysis_data"

if [[ -v IN_SLURM ]]; then
  module load -f midway2 gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
  OUTPUT_ARGS="--output-path=${SLURM_INOUT_DIR}"
  CORES=$SLURM_CPUS_PER_TASK
else
  OUTPUT_ARGS="--output-path=data/stan_analysis_data"
  CORES=8
fi

STAN_THREADS=$((CORES / 4))

echo "Version: $VERSION"
echo "Model: $MODEL"
echo "Iter: $ITER"
echo "Stan threads: $STAN_THREADS"

Rscript --no-save --no-restore --verbose \
  scripts/structural/run-model.R takeup fit \
  --models=${MODEL} \
  --cmdstanr \
  ${OUTPUT_ARGS} \
  --threads=${STAN_THREADS} \
  --outputname=dist_fit${VERSION} \
  --num-mix-groups=1 \
  --chains=4 \
  --iter=${ITER} \
  --sequential
