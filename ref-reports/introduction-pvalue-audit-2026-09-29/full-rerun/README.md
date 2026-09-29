# Completed affected-output rerun — 29 September 2026

**Complete for the affected manuscript outputs; promoted after user approval.** See the [promotion record](../promotion/README.md). The comparisons below use the pre-promotion manuscript snapshots. All **13 affected included tables** were regenerated, and a historical version reproduces the displayed numerical cells of each. The production script produced **26 paired bootstrap model outputs**, plus the level/interval outputs and analytical/simulation appendix tables. All **104 reported sample-size cells** in the paired bootstrap outputs are unchanged.

The reruns retain the existing 500 community-bootstrap draws and, for module attrition, 1,000 imputations. No new inference exercises have been added.

## Main findings

| Result | Before | Corrected | Assessment |
|---|---:|---:|---|
| Continuous-distance pooled take-up interaction p-value | 0.048525 | 0.049492 | Still below 5%, narrowly; both display as 0.049 |
| Continuous-distance pooled observability interaction p-value | 0.002325 | 0.001812 | Remains significant at 1% |
| Binary-distance pooled take-up interaction p-value | 0.082593 | 0.087944 | Still below 10%, not 5% |
| Pooled correct-classification effect in Far communities | 8.50 pp; p=0.003 | 8.58 pp; p=0.003 | Positive and significant |
| Bracelet overall Lee lower bound | 12.2 pp | 12.7 pp | Slightly stronger; these are bound estimates, not confidence limits |
| Preference robustness | No significance-threshold changes | No significance-threshold changes | No new evidence of a distance gradient in private relative valuation |

The main take-up and observability outputs agree exactly with the independent corrected fits from the preceding audit. The main conclusions survive. Some secondary tests cross conventional thresholds; these should be reviewed explicitly rather than described as all results being unchanged.

## Significance and attrition details

[significance-review.md](significance-review.md) lists the **six threshold crossings in the included bootstrap tables**, using unrounded p-values:

- Externality-knowledge pooled average effect: p=0.045 → 0.061 (loses 5%, retains 10%).
- Predicted participation, Bracelet in Close communities: p=0.008 → 0.013 (loses 1%, retains 5%).
- Bracelet binary-distance observability interaction: p=0.012 → 0.009 (gains 1%).
- Conditional accuracy, Calendar in Far communities (decomposition Panel B): p=0.108 → 0.098 (gains 10%).
- Correct classification, Ink distance interaction (Panel C): p=0.045 → 0.058 (loses 5%, retains 10%). The pooled correct-classification interaction remains significant (p=0.016 → 0.015).
- Continuous-distance take-up, Ink interaction: p=0.067 → 0.047 (gains 5%).

The **module-attrition simulation table has one further crossing**: Bracelet's pooled missingness difference changes from −5.83 pp (p=0.111101) to −6.52 pp (p=0.066731), gaining 10% significance. See [module-attrition-contrasts.csv](module-attrition-contrasts.csv), calculated from the saved unrounded simulation output, not rounded table cells.

Module-attrition joint p-values change from 0.431/0.220 to 0.270/0.158; Bracelet–Calendar joint comparisons change from 0.316/0.111 to 0.205/0.119. These remain above 10%. The Close-only Bracelet–Calendar comparison was already significant and stays so (0.041 → 0.039). The Close Bracelet missingness difference becomes −9.1 pp rather than −8.8 pp. The manuscript prose has now been reconciled to these results.

Administrative-sample attrition joint p-values change from 0.081/0.195 to 0.063/0.181. Thus the pooled test remains below 10%, not 5%; this is not a new threshold crossing. Baseline-imbalance observability estimates change modestly without changing their displayed significance categories. The baseline-imbalance take-up table already used the correct original-ID join and is numerically unchanged.

The separate significance-star formatter correction changes annotations at the documented 1%/5%/10% cutoffs; it must not be confused with a change in statistical evidence.

## Which version was in the manuscript?

Ten affected tables reproduce using the 100-draw expected-distance control and old erroneous join. Three later tables reproduce using the 500-draw control and old join: administrative-sample attrition, simulated module attrition, and baseline-imbalance observability. This confirms mixed numerical vintages in the manuscript. [included-table-coverage.csv](included-table-coverage.csv) records these checks.

The two without-controls robustness tables are identical before and after, as expected. SMS, heterogeneity without this control, descriptive recognition/meaning checks, and the structural/policy calculations do not acquire this join correction merely because they appear alongside the affected tables. They were not claimed to be newly re-estimated by this audit.

The broad production script also renders an incentive-check table with an extra drawback column that the manuscript omits. Its five shared columns reproduce the paper. This is an existing presentation difference, unrelated to the join; do not promote that extra column or blindly copy the entire generated directory.

## Review packet

- [manuscript-bootstrap-contrasts.csv](manuscript-bootstrap-contrasts.csv): 276 unrounded contrasts relevant to included bootstrap tables, with the matching manuscript baseline.
- [all-bootstrap-contrasts.csv](all-bootstrap-contrasts.csv): all 824 contrasts, including supporting outputs; separates simulation precision from the join correction using the retained August output where available.
- [significance-review.md](significance-review.md): threshold changes in included bootstrap tables.
- [module-attrition-contrasts.csv](module-attrition-contrasts.csv): exact module-attrition coefficient and p-value comparisons.
- `manuscript/`, `historical/`, `prejoin-500/`, `corrected/`: review copies of tables and aggregate results.
- `table-diffs/`: exact manuscript-versus-corrected LaTeX diffs. These include formatting differences and are **not** an approved patch.
- [completion-checks.json](completion-checks.json), [completed-source-and-result-sha256.json](completed-source-and-result-sha256.json), and [software-versions.txt](software-versions.txt): completion checks, provenance, and package versions.

## Execution and scope

Corrected analysis ran under `build/expected-distance-full-audit-20260929/work/realized`; 100-draw historical reproduction under `build/expected-distance-historical-audit-20260929/work/realized`; additional 500-draw prejoin appendix reconstructions under `build/expected-distance-prejoin-audit-20260929/work/realized`.

The historical launcher first validates a corrected context, then deliberately reconstructs the historical error in memory solely for comparison. Production validation remains enabled. Formulae, sample definitions, distance labels, bootstrap seeds, and inference methods were held fixed. The retained August predicted-participation CSV was also located during this broader review; the earlier audit's statement that it was unavailable is superseded.

The installed dplyr version rejected two obsolete `group_by(..., add=TRUE)` calls. These were changed to `.add=TRUE`, preserving grouping semantics. Endline work then resumed with validated completed outputs reused as checkpoints. Missing preview-directory and historical-launcher completion issues were resolved; successful resumes and section completion markers are retained in the working logs. Both endline resumes and all appendix jobs completed successfully. No manuscript compilation or promotion was performed.

The hash check covers 256 manuscript TeX files. It detected an unrelated concurrent edit to `structural-alternative-models.tex`; this task did not write to Overleaf and left that edit alone. All audited numerical manuscript files remain unchanged.

**Next step: review the numerical replacements and associated prose with the user. The manuscript TODO remains open until those changes are approved and applied.**
