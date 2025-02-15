#!/usr/bin/env bash
#' This file creates various policymaker allocation counterfactuals using the  
#' posterior draws from the structural model. There's a million different arguments 
#' as the script has slowly grown a life of its own over time.

LATEST_VERSION=86 # Version of the model fit to use.
# VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided
VERSION=86

NUM_CORES=16

# Setting arguments
PRED_DISTANCE="" # --pred-distance
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
NUM_POST_DRAWS=200

## What to skip/run
POSTERIOR_MEDIAN="" # --posterior-median or ""
SKIP_PREDICTION=0 # 1
SKIP_OA=0 # 1 or 0
SKIP_PP=0 # 1 or 0
RUN_TARGET_CREATION=0
RUN_ESTIMATION="--run-estimation"

## Welfare function, constraint type, geographic setting
WELFARE_FUNCTION="identity"
CONSTRAINT_TYPE="agg"
COUNTY="full"

## More options...
CUTOFF="" # either no- or empty string - Whether to introduce a hard cutoff to avoid extrapolation
SOLVER="gurobi" # Which solver to use
MANY_POTS="--many-pots" #"--many-pots" - Use all available PoTs or just those in experiment.
SUPPRESS_REP="" # "suppress-rep-" - Remove reputational returns (for some counterfactuals)
CONSTRAINT_TARGET="rep" # Whether to try and hit a target takeup that uses rep returns
STATIC_SIGNAL_PM="" # "--static-signal-pm" - Does the policymaker think the signal is static or dynamic (B&T)
STATIC_SIGNAL_DIST=500 # If static, where to estimate the effect of signals
DEMAND_NAME="" # "static-" - prefix demand filename for whatever reason
DEFAULT_CONSTRAINT_DISTANCE=3500 # Distance cutoff constraint for extrapolation
CONSTRAINT_DISTANCE=${1:-$DEFAULT_CONSTRAINT_DISTANCE} # Get version from command line if provided

## Output Paths
OUTPUT_PATH="optim/data/${MODEL}/${CONSTRAINT_TYPE}-${COUNTY}-many-pots/dist-constraint-${CONSTRAINT_DISTANCE}/" 
PLOT_OUTPUT_PATH="optim/plots/${MODEL}/${CONSTRAINT_TYPE}-${COUNTY}-many-pots" 
DATA_INPUT_NAME="${COUNTY}-many-pots-experiment.rds"

## Fixing \mu or \delta at certain distances
DEFAULT_FIX_PARAM=""
DEFAULT_FIX_PARAM_DISTANCE=""
FIX_PARAM=${2:-${DEFAULT_FIX_PARAM}}
FIX_PARAM_DISTANCE=${3:-${DEFAULT_FIX_PARAM_DISTANCE}}
FIX_PARAM_NAME="${FIX_PARAM}${FIX_PARAM_DISTANCE}"
echo "Constraint distance: ${CONSTRAINT_DISTANCE}"
echo "Fix param: ${FIX_PARAM}"
echo "Fix param distance: ${FIX_PARAM_DISTANCE}"

mkdir -p ${OUTPUT_PATH}
mkdir -p ${PLOT_OUTPUT_PATH}

# Setting up fixing args
FIX_PARAM_ARG=""
if [ $FIX_PARAM == "mu" ]
then 
    FIX_PARAM_ARG="--fix-mu-distance=${FIX_PARAM_DISTANCE}"
fi

if [ $FIX_PARAM == "delta" ]
then
    FIX_PARAM_ARG="--fix-delta-distance=${FIX_PARAM_DISTANCE}"
fi


echo "Fix param arg: ${FIX_PARAM_ARG}"
echo "Fix param name: ${FIX_PARAM_NAME}"

if [ $FIX_PARAM_NAME != "" ]
then 
    FIX_PARAM_NAME="${FIX_PARAM_NAME}-"
fi

set -e



Rscript ./optim/create-distance-data.R \
    --output-name=${DATA_INPUT_NAME} \
    --num-extra-pots=100 \
    --county-subset=${COUNTY} \
    --distance-cutoff=Inf 

# # Create experiment target
if [[ ! -e ${OUTPUT_PATH}/summ-agg-${WELFARE_FUNCTION}-experiment-target-constraint.csv ]]
then
    Rscript ./optim/create-experiment-target.R \
                    --constraint-type=${CONSTRAINT_TYPE} \
                    --welfare-function=${WELFARE_FUNCTION} \
                    --min-cost \
                    --output-path=${OUTPUT_PATH} \
                    --output-basename=summ-${CONSTRAINT_TYPE}-${WELFARE_FUNCTION} \
                    --cutoff-type=cutoff \
                    --data-input-name=full-many-pots-experiment.rds \
                    --posterior-median \
                    --demand-input-path=${OUTPUT_PATH} \
                    --demand-input-filename=pred-demand-dist-fit${VERSION}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv
fi



if [[ ${POSTERIOR_MEDIAN} == "--posterior-median" ]]
then 
    POSTVAR="median"
else
    POSTVAR="post-draws"
fi


run_optim () {

    if [[ ${CUTOFF} == "no-" ]]
    then
        CUTOFF_DIST=10000
    else
        CUTOFF_DIST=3500
    fi

    if [[ ${SUPPRESS_REP} != "" ]]
    then
        SUP_REP_VAR="--suppress-reputation"
    else
        SUP_REP_VAR=""
    fi

    # Gurobi doesn't play nice with renv atm

    if [ $SKIP_PREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    $1 \
                                    $2 \
                                    --output-name=${FIX_PARAM_NAME}${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --to-csv \
                                    --num-post-draws=${NUM_POST_DRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=${CUTOFF_DIST} \
                                    --num-cores=${NUM_CORES} \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --data-input-name=${DATA_INPUT_NAME} \
                                    --output-path=${OUTPUT_PATH} \
                                    --model=${MODEL} \
                                    --single-chain \
                                    ${STATIC_SIGNAL_PM} \
                                    --static-signal-distance=${STATIC_SIGNAL_DIST} \
                                    ${PRED_DISTANCE} \
                                    ${RUN_ESTIMATION} \
                                    ${FIX_PARAM_ARG} \
                                    ${SUP_REP_VAR}
    fi

    if [ $SKIP_OA != 1 ]
    then
    echo "Running optimization"
        Rscript ./optim/optimal_allocation.R  \
                                    ${POSTERIOR_MEDIAN} \
                                    --num-cores=12 \
                                    --min-cost  \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --target-constraint=summ-${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-experiment-target-constraint.csv \
                                    --output-path=${OUTPUT_PATH} \
                                    --output-filename=target-${CONSTRAINT_TARGET}-distconstraint-${CONSTRAINT_DISTANCE}-util-${WELFARE_FUNCTION}-${FIX_PARAM_NAME}${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --input-path=${OUTPUT_PATH}  \
                                    --data-input-path=optim/data \
                                    --data-input-name=${DATA_INPUT_NAME} \
                                    --time-limit=10000 \
                                    --demand-input-filename=pred-demand-dist-fit${VERSION}-${FIX_PARAM_NAME}${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}.csv \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --solver=${SOLVER} \
                                    --distance-constraint=${CONSTRAINT_DISTANCE}

    fi

    if [ $SKIP_PP != 1 ]
    then
        echo "Running postprocessing"
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIOR_MEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --optim-input-a-filename=target-${CONSTRAINT_TARGET}-distconstraint-${CONSTRAINT_DISTANCE}-util-${WELFARE_FUNCTION}-${FIX_PARAM_NAME}${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --data-input-name=${DATA_INPUT_NAME} \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-target-${CONSTRAINT_TARGET}-distconstraint-${CONSTRAINT_DISTANCE}-util-${WELFARE_FUNCTION}-${FIX_PARAM_NAME}${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff \
                                    --pdf-output-path=presentations/optim-takeup-${MODEL}-fig
    fi

}



compare_option () {
    if [[ ${CONSTRAINT_TARGET} == "rep-" ]]
    then
        TMP_REP_VAR_A=""
        TMP_REP_VAR_B=""
    else 
        TMP_REP_VAR_A="suppress-rep-"
        TMP_REP_VAR_B=""
    fi
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIOR_MEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --optim-input-a-filename=target-${CONSTRAINT_TARGET}-${DEMAND_NAME}${TMP_REP_VAR_A}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --optim-input-b-filename=target-${CONSTRAINT_TARGET}-${DEMAND_NAME}${TMP_REP_VAR_B}${CUTOFF}cutoff-b-$3-mu-$4-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --comp-output-basename=target-${CONSTRAINT_TARGET}-${WELFARE_FUNCTION}-${SUPPRESS_REP}${CUTOFF}cutoff-b1-$1-mu1-$2-b2-$3-mu2-$4-${MODEL}-${POSTVAR} \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=target-${CONSTRAINT_TARGET}-${WELFARE_FUNCTION}-${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff \
                                    --pdf-output-path=presentations/takeup-${MODEL}-fig/
}



# run_optim "control" "control" # run control control
# run_optim "control" "bracelet" # counterfactual varying bracelet visibility
run_optim "bracelet" "bracelet" # now bracelet bracelet

 ## now we suppress reputation completely. 
 ## this is because I didn't think of creating a treatment variable with 0 visibility
 ## so we change a global variable woooo

#  SUPPRESS_REP="suppress-rep-"
#  run_optim "control" "control"

#  # # Now we swap to static signalling, fixed at d = 0.5

#  SUPPRESS_REP="" # turn off suppress rep
#  STATIC_SIGNAL_PM="--static-signal-pm" # "--static-signal-pm"
#  STATIC_SIGNAL_DIST=500
#  DEMAND_NAME="static-" 
#  run_optim "control" "bracelet"

# if [[ ${POSTERIOR_MEDIAN} == "--posterior-median" ]]
# then 
#     # Control plots using realised allocation
#     Rscript ./optim/create-presentation-plots.R \
#                                 --constraint-type=agg \
#                                 --welfare-function=${WELFARE_FUNCTION} \
#                                 --min-cost \
#                                 --output-path=${PLOT_OUTPUT_PATH} \
#                                 --output-basename=target-${CONSTRAINT_TARGET}-util-${WELFARE_FUNCTION}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}-${POSTVAR} \
#                                 --cutoff-type=cutoff \
#                                 --data-input-path=optim/data \
#                                 --data-input-name=${DATA_INPUT_NAME} \
#                                 --posterior-median \
#                                 --pdf-output-path=presentations/takeup-${MODEL}-fig \
#                                 --demand-input-path=${OUTPUT_PATH} \
#                                 --demand-input-filename=pred-demand-dist-fit${VERSION}-cutoff-b-control-mu-control-${MODEL}.csv \
#                                 --model=${MODEL} \
#                                 --fit-version=${VERSION} \
#                                 --distance-constraint=${CONSTRAINT_DISTANCE} 

#     # superceded by create-presentation-plots unless you want demand curves
#     # Rscript ./optim/misc-optim-plots.R \
#                                 # --output-path=${OUTPUT_PATH} \
#                                 # --model=${MODEL} \
#                                 # --fit-version=${VERSION} \
#                                 # --distance-constraint=${CONSTRAINT_DISTANCE} \
#                                 # --welfare-function=${WELFARE_FUNCTION} 


#     Rscript ./optim/create-optim-paper-panel.R \
#                                 --output-path=${PLOT_OUTPUT_PATH} \
#                                 --model=${MODEL} \
#                                 --fit-version=${VERSION} \
#                                 --welfare-function=${WELFARE_FUNCTION} \
#                                 --input-path=${OUTPUT_PATH} \
#                                 --distance-constraint=${CONSTRAINT_DISTANCE} 
# else 
#     Rscript ./optim/compare-optim.R \
#         --input-path=${OUTPUT_PATH} \
#         --output-path=${OUTPUT_PATH} \
#         --many-pots \
#         --model=${MODEL} \
#         --welfare-function=${WELFARE_FUNCTION}
# fi
