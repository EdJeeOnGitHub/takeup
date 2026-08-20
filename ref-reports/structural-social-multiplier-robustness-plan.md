# Plan: Structural Social Multiplier Robustness Checks

## Summary

Use the server/raw Stan outputs to test whether the structural robustness specifications recover the paper's social multiplier, not just the ATEs. The implementation should reuse the existing postprocessing path in `scripts/structural/postprocess-roc.R --sm`, then produce a compact paper-style robustness table comparing the multiplier across the structural models already used in the paper.

This should not re-estimate the models unless the needed generated quantities are absent.

## Models to Check

Run the social-multiplier postprocessor for these specifications:

| Check | Fit | Model |
|---|---:|---|
| Main structural model | 105 | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP` |
| No outliers | 105 | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS` |
| Individual distance, community fixed point, individual visibility | 105 | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS` |
| Individual distance, individual fixed point | 105 | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP` |

Also run these as secondary checks if the raw files are available, because they relate to the observability/missingness concern:

| Check | Fit | Model |
|---|---:|---|
| Corrected observability | 106 | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS` |
| SOB observability | 106 | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB` |

If the server uses different fit numbers for the paper tables, use the fit version that generated the current Overleaf robustness table, but keep the exact model names.

## Implementation Steps

1. Confirm raw chain files exist on the server before running anything expensive.

   Expected naming convention:

   ```bash
   data/stan_analysis_data/dist_fit105_<MODEL>-1.csv
   data/stan_analysis_data/dist_fit105_<MODEL>-2.csv
   data/stan_analysis_data/dist_fit105_<MODEL>-3.csv
   data/stan_analysis_data/dist_fit105_<MODEL>-4.csv
   ```

   For each model, inspect the CSV header and verify these generated quantities exist:

   ```text
   cluster_social_multiplier
   cluster_mu_rep
   cluster_mu_rep_deriv
   cluster_delta
   cluster_delta_deriv
   dist_beta_v
   ```

   If these quantities are absent, stop and report that the existing fit cannot answer the social-multiplier robustness question without rerunning generated quantities or the model.

2. Run the existing postprocessor with `--sm`.

   Main example:

   ```bash
   Rscript scripts/structural/postprocess-roc.R 105 \
     --input-path=data/stan_analysis_data \
     --output-path=temp-data/struct-postprocess \
     --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
     --sm 1 2 3 4
   ```

   Repeat for each robustness model.

   Expected outputs:

   ```text
   temp-data/struct-postprocess/rvar_processed_dist_fit105_sm_draws_<MODEL>_1-4.rds
   temp-data/struct-postprocess/rvar_processed_dist_fit105_sm_summ_<MODEL>_1-4.rds
   ```

3. Create one comparison script, for example:

   ```text
   scripts/appendix/struct-social-multiplier-robustness.R
   ```

   The script should:

   - Read each model's `*_sm_summ_*.rds` and `*_sm_draws_*.rds`.
   - Use `sm_rescaled` as the primary social multiplier measure.
   - Before reporting, confirm the sign convention against the existing main-model output from `scripts/structural/render-paper.R` / `temp-data/social-multiplier-decomposition-values.csv`.
   - If the paper reports the multiplier with the opposite sign, transform all models consistently.
   - Restrict the comparison to the experimental support, `roc_distance <= 2500`.
   - Report values at `500m`, `1500m`, and `2500m`.

4. Produce a compact table with this structure:

   ```text
   Specification | Treatment | M(500m) | M(1500m) | M(2500m)
   ```

   Rows should include each robustness specification crossed with:

   ```text
   Control
   Bracelet
   Calendar
   Ink
   ```

   Each cell should report:

   ```text
   posterior median [95% credible interval]
   ```

   The output files should be:

   ```text
   temp-data/struct-social-multiplier-robustness-summary.csv
   presentations/tables/fit105/struct-social-multiplier-robustness-table.tex
   ```

   The `.tex` should be stripped for Overleaf input: no standalone document wrapper and no table notes, matching the style of the existing reduced-form and structural result inputs.

## Validation

The implementing agent should check:

- The main model's reported multiplier matches the existing paper convention after any sign transformation.
- Each model has all four chains represented.
- Each model has all four treatment labels represented.
- Distances `500`, `1500`, and `2500` are present or are selected by nearest available `roc_distance`.
- The robustness models produce the expected number of summary rows by treatment, distance, and variable.
- The new table is not interpreted as an ATE robustness table; it is explicitly a social-multiplier robustness table.

## Assumptions

- The relevant raw Stan outputs are available on the server even where they are incomplete locally.
- The existing generated quantities are sufficient; no model re-estimation is intended.
- The primary object is the rescaled multiplier, `sm_rescaled`, unless comparison to the existing paper output shows that the paper uses a transformed sign convention.
- The paper-facing table should include the four active robustness checks first; observability variants are secondary unless the team decides to include them in the paper.
