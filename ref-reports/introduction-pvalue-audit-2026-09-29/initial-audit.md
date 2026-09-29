# Introduction p-value audit — 29 September 2026

Status: source reconciliation pending. No manuscript changes have been made for this item.

The active manuscript contains three empty p-values. The saved final-distance-classification results provide the values below, but their generated tables differ slightly from the manuscript tables. Do not silently mix results from the two versions.

| Missing quantity | Saved final-classification estimate | Saved p-value |
|---|---:|---:|
| Calendar observability Far–Close interaction | 0.0362217348 | 0.679 |
| Calendar take-up Far–Close interaction | −0.0015772511 | 0.978 |
| Calendar average take-up effect | 0.0270716881 | 0.334 |

## Source checks

- The retained CSVs come from `build/work/realized/temp-data/tidy-rf-tes/`, with corresponding generated tables under `build/realized/presentations/rf-tables/main-specs/`.
- Top-level `temp-data/tidy-rf-tes/` contains different, anchor-band results and is unsuitable for these manuscript claims.
- The observability table agrees with the saved table except for the combined Bracelet coefficient (0.144 in the manuscript versus 0.145 in the generated table) and formatting.
- The take-up table has several small estimate and test differences. For example, Calendar's interaction is −0.003 in the manuscript versus −0.002 in the generated table; the combined Bracelet–Calendar p-value is 0.048 versus 0.046.
- `audited-results.csv` preserves the unrounded estimates and bootstrap standard errors and recomputes two-sided normal p-values using the procedure in `R/reduced-form/functions.R::add_summ_stats`. All 48 displayed p-values match the source CSVs. No p-values were inferred from rounded manuscript cells.
- `sources.json` records source paths and SHA-256 checksums. Both manuscript and generated table snapshots are retained here.

## Surrounding introduction claims

If the saved final-classification results are adopted, reconcile these claims as well:

| Quantity | Current introduction | Saved results |
|---|---|---|
| Control observability Far–Close p-value | 0.056 | 0.023 |
| Ink observability effect in Far p-value | 0.007 | 0.004 |
| Public-signal observability effects in Close p-values | 0.357–0.644 | Bracelet 0.196; Ink 0.341 |
| Bracelet average take-up effect | 7.5 pp | 7.5904 pp (7.6 pp) |
| Bracelet–Calendar average take-up p-value | 0.048 | 0.046 |
| Calendar take-up interaction | −0.3 pp; elsewhere p=0.959 | −0.1577 pp (−0.2 pp); p=0.978 |

The pooled interaction p-values agree: observability 0.009 and take-up 0.083. Bracelet's average take-up p-value also agrees at 0.008. Ink's average take-up p-value is 0.354, consistent with the introduction's two-decimal 0.35.

The predicted-participation claim (12 pp, p=0.002) is not verified: the available older and anchor-band CSVs are not a confirmed source for the currently included table. It should not be updated using those files without further reconciliation.

## Pending decision

The user has been asked whether to reconcile the tables and associated prose to the saved final-classification results or retain the current tables and leave unmatched p-values pending the original result files. This question concerns which numerical version the manuscript should use; no new estimation has been run.
