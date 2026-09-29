# Numerical reconciliation — 29 September 2026

**Status: promoted to the manuscript on 29 September 2026.** Following review and explicit user authorization, synced the 13 affected tables, two additional star-only robustness corrections, and reconciled three prose files. The PDF builds. See the [promotion record](promotion/README.md) for backups, exact diffs, and validation.

## 1. Why the paper and August results differed

The manuscript tables use expected-distance controls calculated with **100 village simulation draws per site-selection draw**. The newer data use **500**. Commit `5340246` (5 May 2026) changed `vill_seeds = 1:100` to `1:500` in `simulate-treatment-assignment/simulate-community-selection.R`. The current simulation file is dated 2 June 2026; the numerical manuscript tables were introduced in Overleaf commit `12bdeca` on 15 April. The May Overleaf edit changed formatting/labels, not estimates.

This is verified, not just inferred from dates: I reconstructed the 100-draw expected distances by selecting village seeds 1–100 from the retained simulation draws, retained the historical join, and reran the three regressions with 500 community bootstrap draws (seeds 1–500). **All 112 numerical cells checked match the manuscript**: 96 coefficients/standard errors across take-up, observability, and predicted participation, plus 16 hypothesis-test p-values in the two main tables. The only normalization needed was `<0.001` versus LaTeX `$<$0.001`. Stars are a separate issue below.

The 100-to-500 change concerns simulation of the expected-distance control, not the number of bootstrap draws or the community distance classification. The 144 final Close/Far labels agree with those in the preexisting covariate CSV.

Rebuilding the current analysis context from current input files reproduces the August covariate, endline, and knowledge-table data exactly. Rerunning the two main regressions reproduces the saved August estimates and standard errors to numerical precision (maximum differences below 3e-15). The faster bootstrap is not the explanation: its predictions were checked against the original regression route for every audited specification.

## 2. A separate preexisting error in the distance-control join

Both the historical and August analyses contain this error. `cluster_id` in the covariate CSV is a sequential rank (1–144); `cluster.id.x` is the original community ID. The expected-distance table contains 150 identified communities, plus an unidentified aggregate row. The old code converted its factor IDs to ranks and joined those to the analysis ranks. The two rank lists have different members, so **all 144 analysis communities receive another community's expected-distance control**.

Concrete example: original community **9**, analysis rank **1**, received community **1**'s expected distance, **1,371.757 m**, rather than its own **1,111.337 m** (see the exact mapping CSV). The same incorrectly joined `mu_d` was then passed from the covariate data to the endline data by original ID.

The local code fix joins by original IDs, enforces a many-to-one join, and validates that each row's control matches the original community ID. Loading an old erroneous cached context now fails with an explicit rebuild instruction. The corrected loader was rebuilt and checked against the independently corrected regression frames. Treatment arms, final distance groups, samples, regression formulas, bootstrap seeds, and inference method were held fixed.

## 3. Significance stars

`prep_tbl()` used `***` for p<0.001, whereas the manuscript table notes specify p<0.01. It also compared already formatted p-value strings. The local fix uses the unrounded two-sided normal p-value from the estimate and bootstrap standard error, with thresholds 0.01/0.05/0.10. A boundary test verifies that rounding across a threshold does not change the star count. This fix changes table annotations, not estimates or standard errors.

## Exact numerical comparison

All effects below are percentage points. The first column reconstructs the tables; introduction p-values are sometimes inconsistent even with those historical tables and should not be treated as their source.

| Outcome / contrast | Paper reconstruction | 500-draw control, old join | Corrected |
|---|---:|---:|---:|
| takeup: bracelet, combined | 7.511 pp; p=0.008 | 7.590 pp; p=0.008 | 7.789 pp; p=0.007 |
| takeup: bracelet - calendar, combined | 4.818 pp; p=0.048 | 4.883 pp; p=0.046 | 5.169 pp; p=0.034 |
| takeup: calendar, far - close | -0.290 pp; p=0.959 | -0.158 pp; p=0.978 | 0.421 pp; p=0.94 |
| takeup: calendar, combined | 2.693 pp; p=0.334 | 2.707 pp; p=0.334 | 2.621 pp; p=0.359 |
| takeup: signal, far - close | 6.772 pp; p=0.083 | 6.775 pp; p=0.083 | 6.756 pp; p=0.088 |
| observability: calendar, far - close | 3.608 pp; p=0.68 | 3.622 pp; p=0.679 | 3.763 pp; p=0.665 |
| observability: signal, far - close | 13.226 pp; p=0.009 | 13.233 pp; p=0.009 | 13.309 pp; p=0.008 |
| predicted: control, far - close | -12.420 pp; p=0.004 | -12.469 pp; p=0.004 | -12.934 pp; p=0.003 |

For predicted participation, the middle column is a fresh reconstruction using the current 500-draw control and old join; it was not available as a saved August output.

All 72 principal contrasts, unrounded estimates, standard errors, and p-values are in [all-result-changes.csv](all-result-changes.csv). It separates the simulation change from the join correction. [historical-reproduction-checks.csv](historical-reproduction-checks.csv) documents the 112 paper-cell matches.

The corrected blank p-values would be **0.665** (Calendar observability interaction), **0.940** (Calendar take-up interaction), and **0.359** (Calendar average take-up effect). These supersede the initially proposed August values of 0.679, 0.978, and 0.334. These values have now been inserted into the paper.

The corrected main results preserve the broad pattern: public signals improve observability more in Far communities, Bracelet increases take-up, Calendar's take-up interaction is near zero, and the pooled take-up interaction remains imprecisely estimated (p=0.088). This statement does not certify every appendix result or the full manuscript.

## Files changed locally

- `R/reduced-form/context.R`: original-ID join and validation of cached/current controls.
- `scripts/checks/test-distance-spec.R`: fix the test's own incorrect join to use original IDs and require nonmissing controls.
- `R/reduced-form/functions.R`: align significance stars with manuscript notes using unrounded p-values.
- `tests/smoke/reduced-form-table-stars.R`: threshold-boundary regression test.
- This audit directory and the manuscript TODO record.

Review-only tables are under [corrected-candidate/rf-tables/main-specs](corrected-candidate/rf-tables/main-specs): the two main tables and predicted participation. They preserve the manuscript layout and change numeric cells and stars only. No candidate prose has been applied. `approved-replacements.json` is the earlier **superseded** August-only plan and must not be applied.

## Full rerun completed

The subsequent [full affected-output rerun](full-rerun/README.md) has now regenerated all 13 affected included tables and reproduced their manuscript baselines. Continuous-distance take-up retains p=0.049, continuous-distance observability retains p=0.002, and the broader packet documents secondary threshold crossings. The outstanding step is user review and consistent manuscript updates, not these reruns.

## What remains before the paper can be called current

1. Review these three distinct changes with the user; agree on the corrected source and consequent prose updates. Do not edit the manuscript until instructed.
2. **Completed:** regenerate and check the other affected included results. See the full rerun packet for continuous-distance, second-order observability, externality knowledge, belief decomposition, Lee bounds, preferences, attrition, and baseline-imbalance checks. Review the associated prose updates before promotion.
3. Check other included table generators for the same significance-star convention. Tables estimated without the affected control are not changed by this join fix, but that alone does not certify their provenance.
4. Stage the agreed replacements consistently and compile the manuscript. **Do not mark the p-value TODO complete yet.**

A code search found no use of this expected-distance covariate in the main structural data builder or Stan model. No structural or policy rerun has been performed or claimed necessary solely on that basis; their numerical validity is outside this reduced-form audit.

## Reproduction and verification

The audit uses the existing two-sided normal test based on 500 Bayesian community-bootstrap draws, matching the table-generating code. This does not add any of the inference exercises that the user excluded from the paper.

Working files, independently reconstructed contexts, and bootstrap draws are under `build/introduction-pvalue-audit-20260929/`. Retained scripts here document the exact reruns; their paths reference that working directory. `reconstruct-expected.py` recreates the 100-draw control; `historical-rerun.R` reproduces the paper; `rerun.R` compares August inputs with the corrected join; `verify-fix.R` checks the production loader and rejects the old cache. The original audit and source snapshots are retained rather than overwritten. `corrected-sources.json` records hashes of current control inputs and changed analysis code.

Checks completed: current input/context agreement; all 144 original-ID mappings checked; corrected loader validated; old cache rejected; distance-specification test passed; star-boundary test passed; original/fast prediction equivalence for all nine regression runs; 112 historical paper cells reproduced. Single-process bootstrap runs were used after a forked attempt stalled; the stalled run was interrupted and its partial outputs were not used.
