SHELL := /usr/bin/env bash
.DEFAULT_GOAL := help

DISTANCE_SPEC ?= assigned
TAKEUP_THREADS ?= 8
BALANCE_SECTIONS ?= all
RI_DRAWS ?= 500
BOOTSTRAP_DRAWS ?= 500
APPENDIX_BOOTSTRAP_DRAWS ?= 1000
APPENDIX_RI_PERMUTATIONS ?= 99999
POLICY_SMOKETEST ?= 1
POLICY_REPLICATES ?= 5
POLICY_SOLVER ?= auto
POLICY_FAST_SOURCE ?= replication/inputs/policy
POLICY_FAST_BUILD ?= build/policy/cluster-weighted
POLICY_MAIN_BUILD ?= build/policy/posterior
POLICY_MODEL_PACKAGE ?= temp-data/policy-model-robustness-complete-20260827-002500
PAPER_FIGURE_BUILD ?= build/paper-figures/$(DISTANCE_SPEC)
POLICY_PAPER_RETURN ?= build/midway-returns/assigned-distance-midway-20260822-160039/policy
POLICY_PAPER_DRAW ?= 337
CANDIDATE_ROOT ?= build/paper-candidate/$(DISTANCE_SPEC)
CANDIDATE_COMPONENT_ROOT ?= build/candidate-components/$(DISTANCE_SPEC)
CANDIDATE_HPC_ROOT ?= build/candidate-hpc/$(DISTANCE_SPEC)
CANDIDATE_HPC_BUNDLE ?=
PAPER_TEX ?= /home/agent/projects/overleaf/overleaf-takeup/ECM ReStud.tex
STRUCTURAL_RENDER_FIT ?= 105
STRUCTURAL_RENDER_INPUT ?= temp-data/struct-postprocess
STRUCTURAL_WORKSPACE ?= build/structural-workspace/main-core-input.RData
export TAKEUP_DISTANCE_SPEC := $(DISTANCE_SPEC)
export TAKEUP_BUILD_SPECS := $(DISTANCE_SPEC)
export TAKEUP_BALANCE_SECTIONS := $(BALANCE_SECTIONS)
export TAKEUP_RI_DRAWS := $(RI_DRAWS)
export TAKEUP_BOOTSTRAP_DRAWS := $(BOOTSTRAP_DRAWS)
export TAKEUP_THREADS
export TAKEUP_STRUCTURAL_WORKSPACE := $(STRUCTURAL_WORKSPACE)

.PHONY: help setup audit analysis-context reduced-form balance balance-tables structural-workspace structural-data structural-fit \
	structural-postprocess compare-distance paper paper-generated paper-full paper-assets \
	structural-render structural-candidate-render structural-paper policy-paper design-paper design-paper-full auxiliary-paper check \
	replication-package paper-audit stan-inventory optimal-policy optimal-policy-legacy \
	policy-fast-predict policy-fast-optimize policy-fast-summarize policy-tables clean-cache \
	structural-theory-figures structural-multiplier-figure policy-main-prepare \
	policy-main-predict policy-main-optimize policy-main-summarize policy-main-render optimal-policy-main \
	balance-ri figures-baseline figures-reduced-form figures-policy-diagnostics \
	figures-structural-diagnostics policy-paper-figures paper-figures-quick check-paper-figures \
	policy-model-table structural-benchmark-recovery-prepare structural-benchmark-recovery-summarize \
	candidate-local candidate-appendix candidate-hpc-export candidate-hpc-import \
	candidate-auxiliary-rf candidate-policy paper-candidate candidate-check \
	distance-comparison-inputs distance-comparison-report distance-comparison

## make policy-model-table    Render the eleven-model optimal-policy robustness table.
policy-model-table:
	Rscript --no-save --no-restore scripts/policy/render-model-scenario-table.R \
		--input-root="$(POLICY_MODEL_PACKAGE)" \
		--output=appendix/structural-robustness/tables/optim-policy-model-scenarios.tex

BENCHMARK_RECOVERY_PATH ?= temp-data/main-core-benchmark-recovery
BENCHMARK_RECOVERY_FITS ?= $(shell find build/structural-fit/assigned -maxdepth 1 -type f -name 'dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-[1-4].csv' 2>/dev/null | sort | paste -sd, -)

## make structural-benchmark-recovery-prepare Generate 50 assigned-design benchmark simulations and the HPC manifest.
structural-benchmark-recovery-prepare:
	Rscript --no-save --no-restore scripts/appendix/generate-main-core-benchmark-recovery.R \
		--workspace=$(STRUCTURAL_WORKSPACE) \
		--fit-csvs="$(BENCHMARK_RECOVERY_FITS)" \
		--output-path="$(BENCHMARK_RECOVERY_PATH)" \
		--replicates=50 --grid-points=21 --seed=20260827

## make structural-benchmark-recovery-summarize Validate all 50 HMC fits and render appendix artifacts.
structural-benchmark-recovery-summarize:
	Rscript --no-save --no-restore scripts/appendix/summarize-main-core-benchmark-recovery.R \
		--input-path="$(BENCHMARK_RECOVERY_PATH)" --expected-replicates=50 \
		--table-path=appendix/structural-robustness/tables/main-core-benchmark-recovery.tex \
		--figure-path=appendix/structural-robustness/figures/main-core-benchmark-likelihood.pdf

help:
	@sed -n 's/^## //p' Makefile

## make setup                  Restore the locked R environment.
setup:
	Rscript --no-save --no-restore -e 'renv::restore(prompt = FALSE)'

## make audit                  Validate and write the Close/Far crosswalk.
audit:
	Rscript --no-save --no-restore scripts/checks/validate-distance-spec.R --spec=$(DISTANCE_SPEC)

## make reduced-form           Build reduced-form outputs for DISTANCE_SPEC.
reduced-form:
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of(c("distance_audit", "build_manifest", "reduced_form")), callr_function = NULL)'

## make analysis-context       Build the shared cleaned context for DISTANCE_SPEC.
analysis-context:
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of("analysis_context"), callr_function = NULL)'

## make balance                Build cached balance sections (BALANCE_SECTIONS=main,orig,...).
balance:
	Rscript --no-save --no-restore -e 'targets::tar_make_future(names = tidyselect::all_of(c("distance_audit", "build_manifest", "balance")), workers = as.integer(Sys.getenv("TAKEUP_THREADS", "1")), callr_function = NULL)'

## make balance-tables         Render tables after all required balance sections exist.
balance-tables:
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of("balance_tables"), callr_function = NULL)'

## make balance-ri             Build only the balance randomization-inference section.
balance-ri:
	TAKEUP_BALANCE_SECTIONS=fit-ri Rscript --no-save --no-restore -e 'targets::tar_make_future(names = tidyselect::all_of(c("distance_audit", "build_manifest", "balance")), workers = as.integer(Sys.getenv("TAKEUP_THREADS", "1")), callr_function = NULL)'

## make structural-workspace   Rebuild minimal structural inputs from current loaders.
structural-workspace: audit
	mkdir -p "$(dir $(STRUCTURAL_WORKSPACE))"
	TAKEUP_DISTANCE_SPEC=realized Rscript --no-save --no-restore \
		scripts/structural/run-model.R takeup fit --data-only \
		--models=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
		--outputname=$(basename $(notdir $(STRUCTURAL_WORKSPACE))) \
		--output-path=$(patsubst %/,%,$(dir $(STRUCTURAL_WORKSPACE))) \
		--num-mix-groups=1 --cmdstanr

## make structural-data        Validate structural inputs for DISTANCE_SPEC.
structural-data: audit structural-workspace
	Rscript --no-save --no-restore scripts/workflow/prepare-structural-distance-data.R --spec=$(DISTANCE_SPEC)

## make structural-fit         Run four latest slim structural chains locally.
structural-fit: structural-data
	bash scripts/workflow/run-structural-fit.sh $(DISTANCE_SPEC)

## make structural-postprocess Rebuild compact GQ from four saved chains.
structural-postprocess: structural-workspace
	Rscript --no-save --no-restore -e 'targets::tar_make(names = tidyselect::all_of("compact_gq"), callr_function = NULL)'

## make structural-render      Render current fit-105 tables/figures from focused summaries for review.
structural-render:
	Rscript --no-save --no-restore scripts/structural/render-paper.R \
		--fit-version=$(STRUCTURAL_RENDER_FIT) \
		--input-path=$(STRUCTURAL_RENDER_INPUT) \
		--output-path=build/structural-paper/$(DISTANCE_SPEC)/fit$(STRUCTURAL_RENDER_FIT)/data \
		--table-output=build/structural-paper/$(DISTANCE_SPEC)/fit$(STRUCTURAL_RENDER_FIT)/tables \
		--figure-output=build/structural-paper/$(DISTANCE_SPEC)/fit$(STRUCTURAL_RENDER_FIT)/figures \
		--write-robustness
	Rscript --no-save --no-restore scripts/workflow/compare-generated-paper-artifacts.R \
		--owner=structural_postprocess \
		--generated-root=build/structural-paper/$(DISTANCE_SPEC)/fit$(STRUCTURAL_RENDER_FIT) \
		--output=build/structural-paper/$(DISTANCE_SPEC)/fit$(STRUCTURAL_RENDER_FIT)/comparison.csv

## make structural-theory-figures Regenerate invariant theoretical Figure 3 panels.
structural-theory-figures:
	Rscript --no-save --no-restore scripts/structural/render-theory-figures.R \
		--output-path=build/structural-paper/invariant/figures

## make structural-multiplier-figure Rebuild and render the focused smooth Figure 4 decomposition.
structural-multiplier-figure: policy-main-prepare
	Rscript --no-save --no-restore scripts/structural/postprocess-sm-curve.R \
		--parameter-csv="$(POLICY_MAIN_BUILD)/policy-posterior-parameters.csv" \
		--output-path=build/structural-paper/invariant/fit105/data/rvar_processed_dist_fit105_sm_summ_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_1-4.rds
	Rscript --no-save --no-restore scripts/structural/render-paper.R \
		--fit-version=105 --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
		--input-path="$(STRUCTURAL_RENDER_INPUT)" \
		--sm-input=build/structural-paper/invariant/fit105/data/rvar_processed_dist_fit105_sm_summ_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_1-4.rds \
		--output-path=build/structural-paper/invariant/fit105/data \
		--table-output=build/structural-paper/invariant/fit105/tables \
		--figure-output=build/structural-paper/invariant/fit105/figures

## make structural-candidate-render Postprocess and render the assigned slim baseline fit.
structural-candidate-render: structural-postprocess
	mkdir -p "$(CANDIDATE_COMPONENT_ROOT)/structural-data"
	Rscript --no-save --no-restore scripts/structural/postprocess-main-core-compact.R \
		--compact-gq-path=build/$(DISTANCE_SPEC)/structural/compact-gq \
		--fit-path=build/structural-fit/$(DISTANCE_SPEC) \
		--output-path="$(CANDIDATE_COMPONENT_ROOT)/structural-data" \
		--fit-version=105 --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
	Rscript --no-save --no-restore scripts/structural/render-paper.R \
		--fit-version=105 --input-path="$(CANDIDATE_COMPONENT_ROOT)/structural-data" \
		--output-path="$(CANDIDATE_COMPONENT_ROOT)/structural-data" \
		--table-output="$(CANDIDATE_COMPONENT_ROOT)/tables" \
		--figure-output="$(CANDIDATE_COMPONENT_ROOT)/figures" --tables-only

## make candidate-appendix       Run isolated local appendix analyses for review.
candidate-appendix: analysis-context
	mkdir -p "$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness"
	Rscript --no-save --no-restore scripts/appendix/campaign-day-original-assignment.R \
		--bootstrap-draws=$(APPENDIX_BOOTSTRAP_DRAWS) \
		--output-dir="$(CANDIDATE_COMPONENT_ROOT)/reviews/campaign-day" \
		--appendix-dir="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness"
	Rscript --no-save --no-restore scripts/appendix/reduced-form-distance-multiplier.R \
		--bootstrap-draws=$(APPENDIX_BOOTSTRAP_DRAWS) \
		--output-dir="$(CANDIDATE_COMPONENT_ROOT)/reviews/reduced-form-distance-multiplier"
	mkdir -p "$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables" \
		"$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/figures"
	cp "$(CANDIDATE_COMPONENT_ROOT)/reviews/reduced-form-distance-multiplier/tables/reduced-form-distance-multiplier.tex" \
		"$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables/reduced-form-distance-multiplier.tex"
	cp "$(CANDIDATE_COMPONENT_ROOT)/reviews/reduced-form-distance-multiplier/figures/pooled-distance-curves.pdf" \
		"$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/figures/reduced-form-pooled-distance-curves.pdf"
	cp "$(CANDIDATE_COMPONENT_ROOT)/reviews/reduced-form-distance-multiplier/figures/derivative-contrast.pdf" \
		"$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/figures/reduced-form-derivative-contrast.pdf"
	Rscript --no-save --no-restore scripts/appendix/randomization-inference.R \
		--permutations=$(APPENDIX_RI_PERMUTATIONS) --cores=$(TAKEUP_THREADS) \
		--output-dir="$(CANDIDATE_COMPONENT_ROOT)/reviews/randomization-inference" \
		--appendix-dir="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness"

## make candidate-auxiliary-rf    Render legacy main-paper RF figures in isolation.
candidate-auxiliary-rf: paper-generated
	Rscript --no-save --no-restore scripts/reduced-form/render-candidate-auxiliary.R \
		--spec=$(DISTANCE_SPEC) --output-root="$(CANDIDATE_COMPONENT_ROOT)"

## make figures-baseline       Render Appendix Figures A1-A2 from the analysis context.
figures-baseline: analysis-context
	Rscript --no-save --no-restore scripts/figures/render-baseline-descriptives.R \
		--context-path=build/$(DISTANCE_SPEC)/context/analysis-context.rds \
		--output-path="$(PAPER_FIGURE_BUILD)/baseline"

## make figures-reduced-form   Render/stage Appendix Figures A8-A10 and E2.
figures-reduced-form: reduced-form balance-ri
	Rscript --no-save --no-restore create-figures/create-gpt-reason-plot.R \
		--output-dir="$(PAPER_FIGURE_BUILD)/reduced-form/figures" \
		--output-basename=gpt-reason-plot --formats=pdf
	Rscript --no-save --no-restore scripts/reduced-form/render-candidate-auxiliary.R \
		--spec=$(DISTANCE_SPEC) --output-root="$(PAPER_FIGURE_BUILD)/reduced-form"

## make figures-policy-diagnostics Render Appendix Figure A11 without optimization.
figures-policy-diagnostics: policy-main-prepare
	Rscript --no-save --no-restore scripts/policy/render-derivative-diagnostic.R \
		--parameter-csv="$(POLICY_MAIN_BUILD)/policy-posterior-parameters.csv" \
		--output-path="$(PAPER_FIGURE_BUILD)/policy" --num-draws=100

## make figures-structural-diagnostics Render Appendix Figure K1 from fit-105 chains.
figures-structural-diagnostics:
	@for chain in 1 2 3 4; do test -f "build/structural-fit/$(DISTANCE_SPEC)/dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-$$chain.csv" || { echo "Missing fit-105 chain $$chain; run make structural-fit DISTANCE_SPEC=$(DISTANCE_SPEC)." >&2; exit 1; }; done
	Rscript --no-save --no-restore scripts/appendix/sigma-u-prior-posterior-plot.R \
		--input-path=build/structural-fit/$(DISTANCE_SPEC) \
		--output-path="$(PAPER_FIGURE_BUILD)/structural"

## make policy-paper-figures   Re-render Figures 5, 6 and A12 from imported 999-draw results.
policy-paper-figures:
	test -f "$(POLICY_PAPER_RETURN)/policy-bootstrap-parameters.csv"
	test -f "$(POLICY_PAPER_RETURN)/allocations/control/replicate-$(shell printf '%04d' $(POLICY_PAPER_DRAW)).rds"
	test -f "$(POLICY_PAPER_RETURN)/allocations/bracelet/replicate-$(shell printf '%04d' $(POLICY_PAPER_DRAW)).rds"
	Rscript --no-save --no-restore scripts/policy/render-paper-figures.R \
		--parameter-csv="$(POLICY_PAPER_RETURN)/policy-bootstrap-parameters.csv" \
		--parameter-draw=median --allocation-draw=$(POLICY_PAPER_DRAW) \
		--policy-path="$(POLICY_PAPER_RETURN)" \
		--output-path=build/policy/cluster-weighted/figures

## make paper-figures-quick    Render all supported quick/local paper figures.
paper-figures-quick: figures-baseline figures-reduced-form figures-policy-diagnostics \
	figures-structural-diagnostics structural-theory-figures structural-multiplier-figure \
	policy-paper-figures
	$(MAKE) check-paper-figures

## make check-paper-figures    Audit manuscript figure classification and generated PDFs.
check-paper-figures:
	Rscript --no-save --no-restore scripts/checks/check-paper-figures.R \
		--paper-tex="$(PAPER_TEX)" --contract=replication/paper-artifact-contract.csv

## make candidate-local          Run all local assigned candidate components.
candidate-local: audit paper-generated structural-fit structural-postprocess structural-candidate-render candidate-appendix candidate-auxiliary-rf

## make candidate-hpc-export     Write assigned-distance external job manifests; submit nothing.
candidate-hpc-export:
	Rscript --no-save --no-restore scripts/workflow/export-candidate-hpc.R \
		--spec=$(DISTANCE_SPEC) --output-root="$(CANDIDATE_HPC_ROOT)/export"

## make candidate-hpc-import     Validate and import an externally completed candidate bundle.
candidate-hpc-import:
	test -n "$(CANDIDATE_HPC_BUNDLE)"
	Rscript --no-save --no-restore scripts/workflow/import-candidate-hpc.R \
		--spec=$(DISTANCE_SPEC) --export-root="$(CANDIDATE_HPC_ROOT)/export" \
		--bundle="$(CANDIDATE_HPC_BUNDLE)" \
		--output-root="$(CANDIDATE_HPC_ROOT)/imported"

## make candidate-policy         Run production sparse policy from assigned imported parameters.
candidate-policy:
	test -f "$(CANDIDATE_HPC_ROOT)/imported/artifacts/policy/policy-bootstrap-parameters.csv"
	$(MAKE) optimal-policy DISTANCE_SPEC=$(DISTANCE_SPEC) \
		POLICY_FAST_SOURCE="$(CANDIDATE_HPC_ROOT)/imported/artifacts/policy" \
		POLICY_FAST_BUILD="$(CANDIDATE_COMPONENT_ROOT)/policy" \
		POLICY_REPLICATES=$(POLICY_REPLICATES)
	mkdir -p "$(CANDIDATE_COMPONENT_ROOT)/reviews/policy-population" \
		"$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables" \
		"$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/figures"
	Rscript --no-save --no-restore scripts/policy/run-population-cost.R \
		--parameter-csv="$(CANDIDATE_HPC_ROOT)/imported/artifacts/policy/policy-bootstrap-parameters.csv" \
		--parameter-type=canonical --analysis-id=exponential-cluster-weights \
		--output-path="$(CANDIDATE_COMPONENT_ROOT)/reviews/policy-population" \
		--work-path="$(CANDIDATE_COMPONENT_ROOT)/policy/population-work" \
		--cores=$(TAKEUP_THREADS) --solver=$(POLICY_SOLVER)
	Rscript --no-save --no-restore scripts/policy/generate-distance-cap-table.R \
		--parameter-csv="$(CANDIDATE_HPC_ROOT)/imported/artifacts/policy/policy-bootstrap-parameters.csv" \
		--parameter-type=canonical \
		--csv-path="$(CANDIDATE_COMPONENT_ROOT)/reviews/policy-population/policy-distance-cap-diagnostics.csv" \
		--table-path="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables/policy-distance-cap-feasibility.tex"
	Rscript --no-save --no-restore scripts/policy/assemble-population-cost.R \
		--input-path="$(CANDIDATE_COMPONENT_ROOT)/reviews/policy-population" \
		--table-path="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables" \
		--figure-path="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/figures"
	Rscript --no-save --no-restore scripts/policy/validate-population-cost.R \
		--input-path="$(CANDIDATE_COMPONENT_ROOT)/reviews/policy-population"
	test -d "$(CANDIDATE_HPC_ROOT)/imported/artifacts/policy/model-robustness"
	Rscript --no-save --no-restore scripts/policy/assemble-model-robustness.R \
		--input-root="$(CANDIDATE_HPC_ROOT)/imported/artifacts/policy/model-robustness" \
		--baseline-csv="$(CANDIDATE_COMPONENT_ROOT)/policy/policy-cluster-bootstrap-summary.csv" \
		--output-csv="$(CANDIDATE_COMPONENT_ROOT)/policy/policy-model-robustness-summary.csv" \
		--table-path="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables/optim-policy-model-robustness.tex" \
		--contrast-table-path="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables/optim-policy-model-robustness-contrasts.tex"
	Rscript --no-save --no-restore scripts/policy/render-model-scenario-table.R \
		--input-root="$(POLICY_MODEL_PACKAGE)" \
		--output="$(CANDIDATE_COMPONENT_ROOT)/appendix/structural-robustness/tables/optim-policy-model-scenarios.tex"

## make paper-candidate          Strictly stage and compile manuscript plus robustness appendix.
paper-candidate:
	Rscript --no-save --no-restore scripts/workflow/stage-candidate-paper.R \
		--spec=$(DISTANCE_SPEC) --stage="$(CANDIDATE_ROOT)" --strict

## make candidate-check          Validate provenance, completeness, isolation, and PDFs.
candidate-check:
	Rscript --no-save --no-restore scripts/checks/check-candidate-paper.R \
		--spec=$(DISTANCE_SPEC) --stage="$(CANDIDATE_ROOT)"

## make compare-distance       Build both definitions and compare artifacts.
compare-distance:
	TAKEUP_BUILD_SPECS=both Rscript --no-save --no-restore -e 'targets::tar_make_future(names = tidyselect::all_of("distance_comparison"), workers = as.integer(Sys.getenv("TAKEUP_THREADS", "1")), callr_function = NULL)'

## make distance-comparison-inputs Build modern realized/assigned review inputs without promotion.
distance-comparison-inputs: compare-distance
	Rscript --no-save --no-restore scripts/figures/render-baseline-descriptives.R \
		--context-path=build/realized/context/analysis-context.rds \
		--output-path=build/paper-figures/realized/baseline
	Rscript --no-save --no-restore scripts/figures/render-baseline-descriptives.R \
		--context-path=build/assigned/context/analysis-context.rds \
		--output-path=build/paper-figures/assigned/baseline
	Rscript --no-save --no-restore create-figures/create-gpt-reason-plot.R \
		--output-dir=build/paper-figures/realized/reduced-form/figures \
		--output-basename=gpt-reason-plot --formats=pdf
	Rscript --no-save --no-restore create-figures/create-gpt-reason-plot.R \
		--output-dir=build/paper-figures/assigned/reduced-form/figures \
		--output-basename=gpt-reason-plot --formats=pdf
	Rscript --no-save --no-restore scripts/reduced-form/render-candidate-auxiliary.R \
		--spec=realized --output-root=build/candidate-components/realized
	Rscript --no-save --no-restore scripts/reduced-form/render-candidate-auxiliary.R \
		--spec=assigned --output-root=build/candidate-components/assigned
	Rscript --no-save --no-restore scripts/reports/build-structural-distance-comparison.R \
		--fit-root=build/structural-fit/assigned --output-root=build/distance-comparison

## make distance-comparison-report Render the paper-order side-by-side audit PDF.
distance-comparison-report:
	Rscript --no-save --no-restore scripts/reports/render-distance-comparison.R \
		--paper-tex="$(PAPER_TEX)" --output-root=build/reports

## make distance-comparison      Build inputs and render the side-by-side audit PDF.
distance-comparison: distance-comparison-inputs distance-comparison-report

## make paper-assets           Validate frozen approved outputs and static assets.
paper-assets:
	Rscript --no-save --no-restore scripts/workflow/manage-paper-artifacts.R

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
	Rscript --no-save --no-restore scripts/workflow/stage-paper.R --spec=$(DISTANCE_SPEC) --strict

## make paper-full             Refresh supported structural and main-policy artifacts, then stage for review.
paper-full: paper-generated structural-postprocess structural-render structural-theory-figures structural-multiplier-figure optimal-policy-main paper-assets
	Rscript --no-save --no-restore scripts/workflow/stage-paper.R --spec=$(DISTANCE_SPEC) --strict --full

## make paper-audit            Report which active paper artifacts are reproducible.
paper-audit:
	Rscript --no-save --no-restore scripts/checks/audit-paper-pipeline-coverage.R

## make stan-inventory         Inventory Stan sources/fits without moving or deleting them.
stan-inventory:
	Rscript --no-save --no-restore scripts/checks/build-stan-artifact-inventory.R

## make policy-tables          Validate synced optimal-policy outputs and render tables.
policy-tables:
	bash archive/code/policy-v1/run_optim_tables.sh

## make policy-fast-predict    Predict sparse feasible-edge demand from compact cluster-weighted parameters.
policy-fast-predict:
	test -f "$(POLICY_FAST_SOURCE)/policy-bootstrap-parameters.csv"
	test -f optim/data/full-many-pots-experiment.rds
	mkdir -p "$(POLICY_FAST_BUILD)"
	Rscript --no-save --no-restore scripts/policy/predict-cluster-bootstrap.R \
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
		Rscript --no-save --no-restore scripts/policy/optimize-cluster-bootstrap.R \
			--input-path="$(POLICY_FAST_BUILD)" \
			--target-csv=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv \
			--scenario-id=$$scenario --num-replicates=$(POLICY_REPLICATES) --solver=$(POLICY_SOLVER); \
	done

## make policy-fast-summarize Render the fast policy table into build/ for review.
policy-fast-summarize: policy-fast-optimize
	Rscript --no-save --no-restore scripts/policy/summarize-cluster-bootstrap.R \
		--input-path="$(POLICY_FAST_BUILD)" \
		--table-path="$(POLICY_FAST_BUILD)/optim-summ-exponential-cluster-weights.tex" \
		--num-replicates=$(POLICY_REPLICATES) --method=exponential

## make policy-main-prepare    Select 200 balanced posterior draws and a median parameter row.
policy-main-prepare:
	Rscript --no-save --no-restore scripts/policy/prepare-baseline-posterior.R \
		--fit-path=build/structural-fit/$(DISTANCE_SPEC) --output-path="$(POLICY_MAIN_BUILD)"

## make policy-main-predict    Predict sparse demand for posterior draws and median parameters.
policy-main-predict: policy-main-prepare
	Rscript --no-save --no-restore scripts/policy/predict-cluster-bootstrap.R \
		--parameter-csv="$(POLICY_MAIN_BUILD)/policy-posterior-parameters.csv" \
		--distance-data=optim/data/full-many-pots-experiment.rds \
		--output-path="$(POLICY_MAIN_BUILD)" --distance-cap=3500 \
		--num-cores=$(TAKEUP_THREADS) --num-replicates=200
	mkdir -p "$(POLICY_MAIN_BUILD)/median"
	Rscript --no-save --no-restore scripts/policy/predict-cluster-bootstrap.R \
		--parameter-csv="$(POLICY_MAIN_BUILD)/policy-posterior-median-parameters.csv" \
		--distance-data=optim/data/full-many-pots-experiment.rds \
		--output-path="$(POLICY_MAIN_BUILD)/median" --distance-cap=3500 \
		--num-cores=1 --num-replicates=1

## make policy-main-optimize   Optimize five scenarios for posterior draws and median parameters.
policy-main-optimize: policy-main-predict
	set -e; for scenario in 1 2 3 4 5; do \
		Rscript --no-save --no-restore scripts/policy/optimize-cluster-bootstrap.R \
			--input-path="$(POLICY_MAIN_BUILD)" \
			--target-csv=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv \
			--scenario-id=$$scenario --num-replicates=200 --solver=$(POLICY_SOLVER); \
		Rscript --no-save --no-restore scripts/policy/optimize-cluster-bootstrap.R \
			--input-path="$(POLICY_MAIN_BUILD)/median" \
			--target-csv=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv \
			--scenario-id=$$scenario --num-replicates=1 --solver=$(POLICY_SOLVER); \
	done

## make policy-main-summarize Render the main posterior policy table.
policy-main-summarize: policy-main-optimize
	Rscript --no-save --no-restore scripts/policy/summarize-cluster-bootstrap.R \
		--input-path="$(POLICY_MAIN_BUILD)" \
		--table-path="$(POLICY_MAIN_BUILD)/optim-summ-table.tex" \
		--num-replicates=200 --method=posterior

## make policy-main-render     Render main-paper demand and allocation-distance figures.
policy-main-render: policy-main-summarize
	Rscript --no-save --no-restore scripts/policy/render-paper-figures.R \
		--median-parameter-csv="$(POLICY_MAIN_BUILD)/policy-posterior-median-parameters.csv" \
		--policy-path="$(POLICY_MAIN_BUILD)/median" \
		--output-path="$(POLICY_MAIN_BUILD)/figures"

## make optimal-policy-main    Run the sparse 200-posterior-draw main-paper policy workflow.
optimal-policy-main: policy-main-render

## make optimal-policy         Run the low-memory sparse policy workflow (5 refits by default).
optimal-policy: policy-fast-summarize

## make optimal-policy-legacy Run the older dense posterior workflow (smoke test by default).
optimal-policy-legacy:
	SMOKETEST=$(POLICY_SMOKETEST) bash archive/code/policy-v1/slurm_run_optim.sh

## make check                  Run build and manuscript dependency checks.
check: audit paper-assets
	Rscript --no-save --no-restore scripts/checks/check-repository-layout.R
	Rscript --no-save --no-restore scripts/checks/test-source-purity.R
	Rscript --no-save --no-restore scripts/checks/test-distance-spec.R
	Rscript --no-save --no-restore scripts/checks/audit-paper-pipeline-coverage.R
	Rscript --no-save --no-restore scripts/checks/build-stan-artifact-inventory.R
	Rscript --no-save --no-restore scripts/checks/check-replication-build.R --spec=$(DISTANCE_SPEC)

## make replication-package    Stage a checksummed replication deposit.
replication-package: check
	Rscript --no-save --no-restore scripts/workflow/stage-replication-package.R

## make clean-cache            Remove only targets metadata and disposable cache.
clean-cache:
	Rscript --no-save --no-restore -e 'targets::tar_destroy(destroy = "meta")'
