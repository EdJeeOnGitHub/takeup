#!/usr/bin/env bash

#SBATCH --partition=broadwl
#SBATCH --job-name=takeup        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=12       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=0-10:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-%j.log
#SBATCH --error=temp/log/takeup-%j.log
#SBATCH --export=IN_SLURM=1
#SBATCH --array=0-2              # robustness batch: 3 tasks

LATEST_VERSION=105
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided
CMDSTAN_ARGS="--cmdstanr"
SLURM_INOUT_DIR="data/stan_analysis_data"
ITER=400

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."
  module load -f midway2 gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0

  OUTPUT_ARGS="--output-path=${SLURM_INOUT_DIR}"
  POSTPROCESS_INOUT_ARGS="--input-path=${SLURM_INOUT_DIR} --output-path=${SLURM_INOUT_DIR}"
  CORES=$SLURM_CPUS_PER_TASK

  echo "Running with ${CORES} cores."
  echo "INOUT ARGS: ${POSTPROCESS_INOUT_ARGS}."
else
  OUTPUT_ARGS="--output-path=data/stan_analysis_data"
  POSTPROCESS_INOUT_ARGS=
  CORES=8
fi

models=(
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS"
)

model=${models[$SLURM_ARRAY_TASK_ID]} # get the model for this task

# stan threads equal to number of cores/4
STAN_THREADS=$((${CORES} / 4))

source quick_postprocess.sh

Rscript --no-save \
  --no-restore \
  --verbose \
  run_takeup.R takeup fit \
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
