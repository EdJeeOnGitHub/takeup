# Observability Selection and Missingness Checks

This note maps the evidence for selection into the knowledge/observability table. The issue is not attrition from the administrative take-up outcome. It is missingness in the endline knowledge table / first-order-beliefs observability module, especially because control households were harder to follow up and may have required more backup respondents.

## 1. Treatment Predicts Knowledge-Table Missingness

What it is: regress an indicator for missing from the knowledge table on treatment assignment, distance assignment, and treatment-by-distance cells.

How it helps: directly tests whether observability-module missingness is different across treatment arms. This is the first diagnostic for whether the observed FOB/observability sample is selected by treatment.

Existing code/table:
- `balance.R`, Knowledge Table Attrition Analysis: `not_in_know_table ~ treatment`.
- Code location: `balance.R:1068`, especially `balance.R:1091`.
- Output: `presentations/tables/attrition-by-treatment.tex`.
- Bootstrap/RF-style version: `scripts/reduced-form/bootstrap.R:1659`.
- Bootstrap outputs: `temp-data/tidy-rf-tes/know-table-attrition-tidy-tes.csv`, table name `know_table_attrition_spec_tbl`.

## 2. Missing vs. Non-Missing Observable Balance

What it is: compare households present in the knowledge table to households missing from it on demographics, distance, baseline knowledge, and implementation/contact measures.

How it helps: tests whether missing households are observably different from observed households. If missing and observed households are similar on pre-treatment and implementation observables, selection into the FOB sample is less concerning.

Existing code/table:
- Code location: `balance.R:1163`.
- Variables listed around `balance.R:1190`.
- Output: `presentations/tables/attrition-covariate-balance.tex`.

## 3. Balance Among Missing Households by Treatment

What it is: restrict to households missing from the knowledge table and test whether their observables differ by treatment arm and distance cell.

How it helps: addresses the sharper selection concern: even if missing households differ from observed households overall, are the missing control households different from missing treatment households? If not, differential treatment-composition among missing cases is less likely to explain treatment effects.

Existing code/table:
- Code location: `balance.R:1299`.
- Output CSV: `temp-data/attrition_treat_comp_df.csv`.
- Output RDS: `temp-data/attrition_treat_joint_tests.rds`.

## 4. Balance Among Observed Households by Treatment

What it is: restrict to households observed in the knowledge table and test whether observables remain balanced by treatment arm and distance cell.

How it helps: the FOB/observability regressions are estimated on observed knowledge-table respondents. This check tests whether that analysis sample remains balanced after missingness.

Existing code/table:
- Code location: `balance.R:1341`.
- Output CSV: `temp-data/non_attrition_treat_comp_df.csv`.
- Output RDS: `temp-data/non_attrition_treat_joint_tests.rds`.

## 5. First-Responder-Only Robustness

What it is: rerun the main FOB/observability treatment-effect table after excluding backup respondents, or equivalently restricting to original sampled/first responders.

How it helps: directly targets the reviewer concern that control households required more backups. If the observability effects are stable when backups are excluded, then the results are not driven by replacing unreachable control respondents with different backup respondents.

Existing code/table:
- Implemented in `scripts/reduced-form/bootstrap.R` in the beliefs regressions section. The code defines `is_backup` from `census_data$endline.backup`, sets `is_first_responder = !is_backup`, writes backup/first-responder sample counts, then reruns the main FOB/observability reduced-form table after `filter(is_first_responder)`.
- Code locations: `scripts/reduced-form/bootstrap.R:1094` for the first-responder flag and sample counts; `scripts/reduced-form/bootstrap.R:1173` for the first-responder-only FOB regression.
- Outputs generated on July 9, 2026 with `devcontainer exec --workspace-folder /home/ed/projects/takeup Rscript scripts/reduced-form/bootstrap.R --beliefs`: `temp-data/tidy-rf-tes/fob-first-responder-sample-counts.csv`, `temp-data/tidy-rf-tes/reducedform-discrete-fob-first-responder-tidy-tes.csv`, `presentations/rf-tables/main-specs/rf_discrete_fob_first_responder_spec_tbl.tex`, and `presentations/rf-tables/main-specs/rf_discrete_fob_first_responder_spec_tbl_weird_order.tex`.
- Main first-responder-only estimates are qualitatively similar to the full FOB table: the signal effect is 0.105 for combined, 0.051 for close, 0.176 for far, and 0.125 for far-minus-close. The bracelet-calendar contrast is 0.104 for combined and 0.127 for far-minus-close.

## 6. Lee Bounds for FOB Observability

What it is: compute Lee-style upper and lower bounds for the FOB/observability estimates using selection rates into the FOB sample.

How it helps: provides a worst-case selection adjustment for missing observability outcomes. This is useful even if missingness is only modestly differential, because it bounds how much the FOB estimates could move under monotone selection.

Existing code/table:
- Code location: `scripts/reduced-form/bootstrap.R:1178`.
- Outputs: `temp-data/tidy-rf-tes/fob-lee-upper-tidy-tes.csv`, `temp-data/tidy-rf-tes/fob-lee-lower-tidy-tes.csv`.
- Table name: `fob_lee_bounds_tbl`.
- Note: the code currently has `stop()` after saving the Lee-bounds table at `scripts/reduced-form/bootstrap.R:1302`, so running the full beliefs section will halt there.

## 7. Calendar Active-Control Checks

What it is: inspect whether calendar has similar missingness/backup patterns and whether calendar-control observability comparisons survive the missingness checks.

How it helps: calendar is useful as an active-control argument. If calendar has comparable follow-up/missingness issues but does not generate the same observability pattern as bracelet, then differential follow-up alone is less likely to explain the substantive observability results.

Existing code/table:
- Treatment-specific attrition estimates are included in `balance.R:1091` and the bootstrapped missingness outputs in `scripts/reduced-form/bootstrap.R:1659`.
- Main calendar/control observability estimates are in the beliefs/FOB reduced-form section of `scripts/reduced-form/bootstrap.R`.
- Needed addition: explicitly pull out and narrate calendar missingness and calendar-control FOB comparisons together.

## 8. Structural Robustness to Missing Beliefs

What it is: rerun or perturb the structural model using alternative assumptions about missing observability/beliefs inputs, such as imputation/sampling over missing observations or alternative observability measures.

How it helps: this is the downstream robustness for the structural estimates. The reduced-form checks diagnose selection; structural robustness shows that the multiplier/decomposition is not fragile to plausible missingness adjustments.

Existing code/table:
- Not identified as already implemented.
- Related TODO in `ref-reports/edited-feedback-refine.md`: sample over missing data for people missing beliefs in the structural model, ideally without a full selection model unless the RF selection checks show one is needed.
