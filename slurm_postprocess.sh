#!/usr/bin/env bash

#SBATCH --partition=bigmem2
#SBATCH --job-name=quick-takeup        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1              # total number of tasks across all nodes
#SBATCH --cpus-per-task=1      # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=50G         # memory per cpu-core (4G is default)
#SBATCH --time=0-10:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-%j.log
#SBATCH --error=temp/log/takeup-%j.log
#SBATCH --array=0-4 # Adjust this to the number of models - 1
#SBATCH --export=IN_SLURM=1

LATEST_VERSION=96
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided
SLURM_INOUT_DIR="data/stan_analysis_data"

models=(
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_SD_WTP_VAL"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_MU_WTP_VAL"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_WTP_SUBMODEL"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_BELIEFS_SUBMODEL"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_SUBMODELS"
)

prior_args=(
  "--prior"
  ""
)


echo "Version: $VERSION"

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."
	

  module load -f midway2 gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
  module load -f R/4.2.0
	IN_ARG="--input-path=${SLURM_INOUT_DIR}"
	OUT_ARG="--output-path=${SLURM_INOUT_DIR}"

  echo "Running with ${CORES} cores."
  echo "INOUT ARGS: ${IN_ARG} ${OUT_ARG}."
else
	IN_ARG="--input-path=data/stan_analysis_data"
	OUT_ARG="--output-path=temp-data"
  CORES=8
fi


model=${models[$SLURM_ARRAY_TASK_ID]} # get the model for this task

Rscript --no-save \
        --no-restore \
        --verbose \
        quick_ate_postprocess.R \
        ${VERSION} \
        --model=${model} \
        1 2 3 4 > temp/log/struct-postprocess_${model}_${VERSION}.txt 2>&1 

	  # Within SLURM tasks
#	  srun --export=all --exclusive --ntasks=1 bash -c \
#	    "source quick_postprocess.sh && postprocess_model \
#	      quick_ate_postprocess.R \
#	      ${VERSION} \
#	      ${model} \
#	      ${IN_ARG} \
#	      ${OUT_ARG} \
#	      ${prior_arg}" &
#	  srun --export=all --exclusive --ntasks=1 bash -c \
#	    "source quick_postprocess.sh && postprocess_model \
#	      quick_submodel_postprocess.R \
#	      ${VERSION} \
#	      ${model} \
#	      ${IN_ARG} \
#	      ${OUT_ARG} \
#	      ${prior_arg}" &
#	  srun --export=all --exclusive --ntasks=1 bash -c \
#	    "source quick_postprocess.sh && postprocess_model \
#	      quick_roc_postprocess.R \
#	      ${VERSION} \
#	      ${model} \
#	      ${IN_ARG} \
#	      ${OUT_ARG} \
#	      ${prior_arg} \
#	      --cluster-roc \
#	      --cluster-takeup-prop \
#	      --cluster-rep-return-dist \
#	      --sm" &
