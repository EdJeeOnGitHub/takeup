SHELL := /usr/bin/env bash
.DEFAULT_GOAL := help

DISTANCE_SPEC ?= assigned
TAKEUP_THREADS ?= 8
BALANCE_SECTIONS ?= all
RI_DRAWS ?= 500
BOOTSTRAP_DRAWS ?= 500
POLICY_SMOKETEST ?= 1
export TAKEUP_DISTANCE_SPEC := $(DISTANCE_SPEC)
export TAKEUP_BUILD_SPECS := $(DISTANCE_SPEC)
export TAKEUP_BALANCE_SECTIONS := $(BALANCE_SECTIONS)
export TAKEUP_RI_DRAWS := $(RI_DRAWS)
export TAKEUP_BOOTSTRAP_DRAWS := $(BOOTSTRAP_DRAWS)
export TAKEUP_THREADS

.PHONY: help setup audit reduced-form balance balance-tables structural-data structural-fit \
	structural-postprocess compare-distance paper paper-generated paper-full paper-assets \
	structural-paper policy-paper design-paper design-paper-full auxiliary-paper check \
	replication-package paper-audit stan-inventory optimal-policy policy-tables clean-cache

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

## make paper-full             Refresh expensive structural/policy intermediates, then stage for review.
paper-full: paper-generated structural-postprocess optimal-policy policy-tables paper-assets
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

## make optimal-policy         Run policy workflow (smoke test by default; POLICY_SMOKETEST=0 for production).
optimal-policy:
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
