#!/usr/bin/env bash

run_optimal_incentives () {
    f_VERSION=$1
    f_MODEL=$2
    f_B=$3
    f_MU=$4
    f_LAMBDA=$5
    f_EXTERNALITY=$6
    f_POSTERIOR=$7

    echo "RUNNING VERSION: $f_VERSION"
    echo "RUNNING MODEL: $f_MODEL"
    echo "RUNNING B: $f_B"
    echo "RUNNING MU: $f_MU"
    echo "RUNNING LAMDBA: $f_LAMBDA"
    echo "RUNNING EXTERNALITY: $f_EXTERNALITY"
    if [ "$f_POSTERIOR" == "--posterior" ]; then
        echo "RUNNING FULL POSTERIOR"
    fi

    Rscript --no-save \
            --no-restore \
            --verbose \
            scratch/optimal-incentives.R \
            $f_VERSION \
            $f_B \
            $f_MU \
	    --model=$f_MODEL \
            --num-post-draws=200 \
            $f_POSTERIOR \
            --robust-externality \
            --robust-lambda \
            --output-name="optimal-incentives-b-$f_B-mu-$f_MU-lambda-$f_LAMBDA-externality-$f_EXTERNALITY-$f_MODEL-fitversion$f_VERSION" \
            > temp/log/optimal-incentives-b-$f_B-mu-$f_MU-lambda-$f_LAMBDA-externality-$f_EXTERNALITY-$f_MODEL-fitversion$f_VERSION.log 2>&1 &
   wait
}

# Run full posterior across all treatments once
run_optimal_incentives 86 "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP" "control" "control" 0 0 "--posterior"
# run just posterior mean part of script for rest
treatments=("ink" "calendar" "bracelet")
for treatment in "${treatments[@]}"
do
   run_optimal_incentives 86 "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP" $treatment $treatment 0 0
done