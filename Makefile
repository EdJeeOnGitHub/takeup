SHELL := /usr/bin/env bash
.DEFAULT_GOAL := help

DISTANCE_SPEC ?= assigned
TAKEUP_THREADS ?= 8
BALANCE_SECTIONS ?= all
RI_DRAWS ?= 500
BOOTSTRAP_DRAWS ?= 500
POLICY_SMOKETEST ?= 1
POLICY_REPLICATES ?= 5
POLICY_SOLVER ?= auto
POLICY_FAST_SOURCE ?= optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots-exponential-cluster-weights
POLICY_FAST_BUILD ?= build/policy/cluster-weighted
STRUCTURAL_RENDER_FIT ?= 105
STRUCTURAL_RENDER_INPUT ?= temp-data/struct-postprocess
export TAKEUP_DISTANCE_SPEC := $(DISTANCE_SPEC)
export TAKEUP_BUILD_SPECS := $(DISTANCE_SPEC)
export TAKEUP_BALANCE_SECTIONS := $(BALANCE_SECTIONS)
export TAKEUP_RI_DRAWS := $(RI_DRAWS)
export TAKEUP_BOOTSTRAP_DRAWS := $(BOOTSTRAP_DRAWS)
export TAKEUP_THREADS

.PHONY: help setup audit reduced-form balance balance-tables structural-data structural-fit \
	structural-postprocess compare-distance paper paper-generated paper-full paper-assets \
	structural-render structural-paper policy-paper design-paper design-paper-full auxiliary-paper check \
	replication-package paper-audit stan-inventory optimal-policy optimal-policy-legacy \
	policy-fast-predict policy-fast-optimize policy-fast-summarize policy-tables clean-cache

help:
	@sed -n 's/^## //p' Makefile

## make setup                  Restore the locked R environment.
setup:
	Rscript --no-save --no-restore -e 'renv::restore(prompt = FALSE)'

## make audit                  Validate and write the Close/Far crosswalk.
audit:
	Rscript --no-save --no-restore scripts/validate-distance-spec.R --spec=$(DISTANCE_SPEC)

## make reduced-form           Build reduced-form outputs for DISTANCE_SPEC.
reduced-form:
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of(c("distance_audit", "build_manifest", "reduced_form")), callr_function = NULL)'

## make balance                Build cached balance sections (BALANCE_SECTIONS=main,orig,...).
balance:
	Rscript --no-save --no-restore -e 'targets::tar_make_future(names = tidyselect::all_of(c("distance_audit", "build_manifest", "balance")), workers = as.integer(Sys.getenv("TAKEUP_THREADS", "1")), callr_function = NULL)'

## make balance-tables         Render tables after all required balance sections exist.
balance-tables:
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of("balance_tables"), callr_function = NULL)'

## make structural-data        Validate structural inputs for DISTANCE_SPEC.
structural-data: audit
	Rscript --no-save --no-restore scripts/prepare-structural-distance-data.R --spec=$(DISTANCE_SPEC)

## make structural-fit         Run four latest slim structural chains locally.
structural-fit: structural-data
	bash scripts/run-structural-fit.sh $(DISTANCE_SPEC)

## make structural-postprocess Rebuild compact GQ from four saved chains.
structural-postprocess:
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of("compact_gq"), callr_function = NULL)'

## make structural-render      Render current fit-105 tables/figures from focused summaries for review.
structural-render:
	Rscript --no-save --no-restore structural_tables.R \
		--fit-version=$(STRUCTURAL_RENDER_FIT) \
		--input-path=$(STRUCTURAL_RENDER_INPUT) \
		--output-path=build/structural-paper/fit$(STRUCTURAL_RENDER_FIT)/data \
		--table-output=build/structural-paper/fit$(STRUCTURAL_RENDER_FIT)/tables \
		--figure-output=build/structural-paper/fit$(STRUCTURAL_RENDER_FIT)/figures \
		--write-robustness
	Rscript --no-save --no-restore scripts/compare-generated-paper-artifacts.R \
		--owner=structural_postprocess \
		--generated-root=build/structural-paper/fit$(STRUCTURAL_RENDER_FIT) \
		--output=build/structural-paper/fit$(STRUCTURAL_RENDER_FIT)/comparison.csv

## make compare-distance       Build both definitions and compare artifacts.
compare-distance:
	TAKEUP_BUILD_SPECS=both Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of("distance_comparison"), callr_function = NULL)'

## make paper-assets           Validate frozen approved outputs and static assets.
paper-assets:
	Rscript --no-save --no-restore scripts/manage-paper-artifacts.R

## make structural-paper      Validate approved structural render outputs.
structural-paper: paper-assets

## make policy-paper          Validate approved optimal-policy outputs.
policy-paper: paper-assets

## make design-paper          Validate approved design simulation/map outputs.
design-paper: paper-assets

## make auxiliary-paper       Validate approved legacy reduced-form/design outputs.
auxiliary-paper: paper-assets

## make design-paper-full     Run the design simulation source workflow, then review before promotion.
design-paper-full:
	Rscript --no-save --no-restore simulate-treatment-assignment/simulate-community-selection.R
	Rscript --no-save --no-restore -e 'dir.create("build/design-paper", recursive = TRUE, showWarnings = FALSE); rmarkdown::render("takeup_pap.Rmd", output_dir = "build/design-paper")'

## make paper-generated       Regenerate the ordinary reduced-form and balance outputs.
paper-generated:
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of(c("distance_audit", "build_manifest", "reduced_form")), callr_function = NULL)'
	Rscript --no-save --no-restore -e 'targets::tar_make_future(names = tidyselect::all_of(c("distance_audit", "build_manifest", "balance")), workers = as.integer(Sys.getenv("TAKEUP_THREADS", "1")), callr_function = NULL)'
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of("balance_tables"), callr_function = NULL)'

## make paper                  Build/stage every paper dependency without replacing approved outputs.
paper: paper-generated paper-assets structural-paper policy-paper design-paper auxiliary-paper
	Rscript --no-save --no-restore scripts/stage-paper.R --spec=$(DISTANCE_SPEC) --strict

## make paper-full             Refresh structural renders and fast policy intermediates, then stage for review.
paper-full: paper-generated structural-postprocess structural-render optimal-policy paper-assets
	Rscript --no-save --no-restore scripts/stage-paper.R --spec=$(DISTANCE_SPEC) --strict

## make paper-audit            Report which active paper artifacts are reproducible.
paper-audit:
	Rscript --no-save --no-restore scripts/audit-paper-pipeline-coverage.R

## make stan-inventory         Inventory Stan sources/fits without moving or deleting them.
stan-inventory:
	Rscript --no-save --no-restore scripts/build-stan-artifact-inventory.R

## make policy-tables          Validate synced optimal-policy outputs and render tables.
policy-tables:
	bash run_optim_tables.sh

## make policy-fast-predict    Predict sparse feasible-edge demand from compact cluster-weighted parameters.
policy-fast-predict:
	test -f "$(POLICY_FAST_SOURCE)/policy-bootstrap-parameters.csv"
	test -f optim/data/full-many-pots-experiment.rds
	mkdir -p "$(POLICY_FAST_BUILD)"
	Rscript --no-save --no-restore optim/predict-policy-cluster-bootstrap.R \
		--parameter-csv="$(POLICY_FAST_SOURCE)/policy-bootstrap-parameters.csv" \
		--distance-data=optim/data/full-many-pots-experiment.rds \
		--output-path="$(POLICY_FAST_BUILD)" --distance-cap=3500 \
		--num-cores=$(TAKEUP_THREADS) --num-replicates=$(POLICY_REPLICATES)

## make policy-fast-optimize   Optimize five policy scenarios from sparse predicted demand.
policy-fast-optimize: policy-fast-predict
	test -f optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv
	@if ! command -v glpsol >/dev/null 2>&1 && ! command -v gurobi_cl >/dev/null 2>&1 && [ ! -x /opt/gurobi1000/linux64/bin/gurobi_cl ]; then \
		echo "No MILP solver found. Install glpk-utils (glpsol) or provide gurobi_cl." >&2; \
		exit 1; \
	fi
	set -e; for scenario in 1 2 3 4 5; do \
		Rscript --no-save --no-restore optim/optimize-policy-cluster-bootstrap.R \
			--input-path="$(POLICY_FAST_BUILD)" \
			--target-csv=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv \
			--scenario-id=$$scenario --num-replicates=$(POLICY_REPLICATES) --solver=$(POLICY_SOLVER); \
	done

## make policy-fast-summarize Render the fast policy table into build/ for review.
policy-fast-summarize: policy-fast-optimize
	Rscript --no-save --no-restore optim/summarize-policy-cluster-bootstrap.R \
		--input-path="$(POLICY_FAST_BUILD)" \
		--table-path="$(POLICY_FAST_BUILD)/optim-summ-exponential-cluster-weights.tex" \
		--num-replicates=$(POLICY_REPLICATES) --method=exponential

## make optimal-policy         Run the low-memory sparse policy workflow (5 refits by default).
optimal-policy: policy-fast-summarize

## make optimal-policy-legacy Run the older dense posterior workflow (smoke test by default).
optimal-policy-legacy:
	SMOKETEST=$(POLICY_SMOKETEST) bash slurm_run_optim.sh

## make check                  Run build and manuscript dependency checks.
check: audit paper-assets
	Rscript --no-save --no-restore scripts/test-distance-spec.R
	Rscript --no-save --no-restore scripts/audit-paper-pipeline-coverage.R
	Rscript --no-save --no-restore scripts/build-stan-artifact-inventory.R
	Rscript --no-save --no-restore scripts/check-replication-build.R --spec=$(DISTANCE_SPEC)

## make replication-package    Stage a checksummed replication deposit.
replication-package: check
	Rscript --no-save --no-restore scripts/stage-replication-package.R

## make clean-cache            Remove only targets metadata and disposable cache.
clean-cache:
	Rscript --no-save --no-restore -e 'targets::tar_destroy(destroy = "meta")'
