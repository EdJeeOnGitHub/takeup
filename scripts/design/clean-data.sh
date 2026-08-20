#!/bin/bash

set -e

Rscript --no-save \
        --no-restore \
        --verbose \
        ./rct-design-fieldwork/clean-baseline-data.R

Rscript --no-save \
        --no-restore \
        --verbose \
        ./rct-design-fieldwork/clean-endline-data.R