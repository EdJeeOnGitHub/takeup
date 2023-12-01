#!/bin/bash
#PBS -N "takeup-pp" 
#PBS -j oe
#PBS -V
#PBS -e /home/edjee/projects/takeup/logs/error.txt
#PBS -o /home/edjee/projects/takeup/logs/output.txt
#PBS -l procs=16,mem=80gb
cd $PBS_O_WORKDIR




fit_models () {

AKARING_OUTPUT_PATH="/share/akaringlab/takeup-output/"

LATEST_VERSION=96
CMDSTAN_ARGS="--cmdstanr --include-paths=stan_models"
VERSION=${1:-$LATEST_VERSION}
#OUTPUT_ARGS="--output-path=$PBS_O_WORKDIR/data/stan_analysis_data"
OUTPUT_ARGS="--output-path=${AKARING_OUTPUT_PATH}/data/stan_analysis_data"
POSTPROCESS_INOUT_ARGS="--input-path=${AKARING_OUTPUT_PATH}/data/stan_analysis_data --output-path=${AKARING_OUTPUT_PATH}/temp-data"
STAN_THREADS=4
ITER=800

models=(
#	"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP" 
#	"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_SD_WTP_VAL"
#	"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_MU_WTP_VAL"
	"REDUCED_FORM_NO_RESTRICT_DIST_CTS"
)

	for model in "${models[@]}"
	do
		Rscript --no-save \
			--no-restore \
			--verbose \
			run_takeup.R takeup fit \
			--models=${model} \
			${CMDSTAN_ARGS} \
			${OUTPUT_ARGS} \
			--threads=${STAN_THREADS} \
			--outputname=dist_fit${VERSION} \
			--num-mix-groups=1 \
			--iter=${ITER} \
			--sequential > logs/output-${model}-fit${VERSION}.txt 2>&1

		Rscript --no-save \
			--no-restore \
			--verbose \
			run_takeup.R takeup prior \
			--models=${model} \
			${CMDSTAN_ARGS} \
			${OUTPUT_ARGS} \
			--threads=${STAN_THREADS} \
			--outputname=dist_prior${VERSION} \
			--num-mix-groups=1 \
			--iter=${ITER} \
			--sequential > logs/output-${model}-prior${VERSION}.txt 2>&1
	done
	

#	Rscript --no-save \
#		--no-restore \
#		--verbose \
#		quick_roc_postprocess.R \
#		${VERSION} \
#		${POSTPROCESS_INOUT_ARGS} \
#		--model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB \
#		--cluster-roc \
#		--cluster-takeup-prop \
#		--cluster-rep-return-dist \
#		--sm \
#		1 2 3 4 > logs/postprocess-output${VERSION}.txt 2>&1 
#
#	Rscript --no-save \
#		--no-restore \
#		--verbose \
#		quick_roc_postprocess.R \
#		${VERSION} \
#		${POSTPROCESS_INOUT_ARGS} \
#		--model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB \
#		--prior \
#		--cluster-roc \
#		--cluster-takeup-prop \
#		--cluster-rep-return-dist \
#		--sm \
#		1 2 3 4 > logs/postprocess-output${VERSION}-prior.txt 2>&1 

	Rscript --no-save \
		--no-restore \
		--verbose \
		quick_ate_postprocess.R \
		${VERSION} \
		${POSTPROCESS_INOUT_ARGS} \
		--model=REDUCED_FORM_NO_RESTRICT_DIST_CTS \
		1 2 3 4 > logs/postprocess-output${VERSION}-model 2>&1 

	Rscript --no-save \
		--no-restore \
		--verbose \
		quick_ate_postprocess.R \
		${VERSION} \
		${POSTPROCESS_INOUT_ARGS} \
		--model=REDUCED_FORM_NO_RESTRICT_DIST_CTS \
		--prior \
		1 2 3 4 > logs/postprocess-output${VERSION}-prior.txt 2>&1 


}

export -f fit_models

# call mpi
mpirun -n 1 -machinefile $PBS_NODEFILE  bash -c 'fit_models'
	
