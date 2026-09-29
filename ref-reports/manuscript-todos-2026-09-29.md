# Manuscript TODOs — 29 September 2026

Review of the active manuscript in `/home/ed/projects/overleaf/overleaf-takeup/`, checked against repository TODOs and September review notes. The main outstanding work is reconciliation and finalization; several older TODO lists are stale.

Scope: sources, included tables, and the existing build log were inspected. Estimates were not rerun and a fresh compilation was not performed. Items below distinguish corrections from decisions about manuscript scope.

## Highest-priority substantive work

- [x] **Align distance terminology with the revised randomization discussion.** Completed 29 September 2026. The revised design discussion remains the reference point. Applied the approved contribution wording, structural distance definition, structural table notes, and PAP analysis description and distance clarification. After review, retained the existing abstract, introduction, and design-overview wording. No anchor-band comparisons are included in the active manuscript. The user closed this manuscript item; the broader code cleanup was not performed. Sources: [design](</home/ed/projects/overleaf/overleaf-takeup/experimental-design.tex:31>), [PAP discussion](</home/ed/projects/overleaf/overleaf-takeup/online-appendix.tex:968>).

- [x] **Use one definition of the social multiplier.** Completed 29 September 2026. Replaced “absent social image incentives” with “holding social image returns fixed” in the [introduction](</home/ed/projects/overleaf/overleaf-takeup/ECM ReStud.tex:271>) and [theory introduction](</home/ed/projects/overleaf/overleaf-takeup/theory.tex:8>). Both are direct replacements; the surrounding text and structural definition already use the intended benchmark.

  **Rationale:** The denominator is the local travel-cost response at the existing equilibrium, holding the social-image return fixed at its prevailing value. Removing social image altogether would generally change the equilibrium cutoff and the density at which the response is evaluated. Bénabou and Tirole's *Laws and Norms* (2011), Section 1.2, equation (6), defines the multiplier as `1 / [1 + μ Δ′(v*)]`, multiplying the density at the equilibrium cutoff. Our expression, `[1 − μ′(d) Δ(w*)] / [1 + μ(d) Δ′(w*)]`, extends this construction to distance-dependent observability and reduces to theirs when `μ′(d) = 0`. The 2025 version, equation (6), retains the local-derivative construction with a different normalization. References: [2011 version, pp. 6–7](https://www.tse-fr.eu/sites/default/files/medias/doc/by/tirole/laws_and_norms_oct3.pdf), [2025 version, p. 8](https://www.tse-fr.eu/sites/default/files/TSE/documents/doc/by/tirole/laws_and_norms_january_21_2025.pdf). Policy scenarios explicitly setting social-image returns to zero are separate counterfactuals and retain their existing labels.

- [x] **Review headline claims against the robustness results.** Reviewed 29 September 2026; no manuscript changes needed. The abstract and introduction establish the structural context, the robustness discussion acknowledges weaker evidence under alternative information and observability assumptions, and the conclusion explains the role of structure in its limitations. After reviewing the specific passages, the user chose to retain the current wording. Sources: [robustness discussion](</home/ed/projects/overleaf/overleaf-takeup/structural-model-rewrite.tex:236>), [conclusion](</home/ed/projects/overleaf/overleaf-takeup/conclusion.tex:8>).

- [x] **Resolve the inference presentation — do not include.** Decision recorded 29 September 2026: the user explicitly excludes main-interaction randomization inference (including multiple-testing adjustments) and structural cluster-weighted uncertainty checks from the manuscript. Retain the underlying analyses as supporting records; do not add them to the paper or reopen their integration as an outstanding TODO. Sources: [RI results](../appendix/structural-robustness/tables/randomization-inference.tex), [completed cluster-weight work](remaining-todos.md).

- [x] **Confirm the scope of policy robustness reporting — intentional omission.** Decision recorded 29 September 2026: the user confirms that the alternative-model policy table is intentionally omitted from the manuscript. Retain the current distance-cap and resource-cost checks; do not treat integration of the alternative-model policy table as outstanding work. Sources: [policy appendix](</home/ed/projects/overleaf/overleaf-takeup/online-appendix.tex:1204>), [earlier promotion plan](../doc/policy-population-manuscript-swap-plan-2026-09-07.md).

- [x] **Qualify the policy forecasting interpretation.** Completed 29 September 2026. The approved wording preserves the planner forecasting-error interpretation and explicitly specifies a planner who uses the estimated private payoffs but ignores social image returns. It does not describe a separately fitted model without social image. Source: [policy interpretation](</home/ed/projects/overleaf/overleaf-takeup/optimal-policy-rewrite.tex:127>).

Earlier scope decisions explicitly excluded pooling, short-cap, and policy-bootstrap exercises. These are not mandatory unfinished manuscript work.

## Concrete corrections and finalization

- [x] **Fill three empty introduction p-values and reconcile numerical results.** Completed 29 September 2026 after the user approved promotion of the full rerun. Synced 13 corrected tables, corrected significance annotations in two numerically unchanged robustness tables, and reconciled statistics and inference language in the introduction, experimental-design section, and online appendix. The three p-values are 0.665, 0.940, and 0.359. Verified every promoted table data cell against the audited outputs and rebuilt the PDF. The existing broken PAP reference and oversized floats remain separate tasks below. The audit separates the 100-to-500 simulation change, original-ID join correction, and significance-star correction. See the [reconciliation audit](introduction-pvalue-audit-2026-09-29/README.md), [full rerun](introduction-pvalue-audit-2026-09-29/full-rerun/README.md), and [promotion record](introduction-pvalue-audit-2026-09-29/promotion/README.md).

- [x] **Complete the alternative-model appendix review.** Resolved by the user on 29 September 2026.

- [x] **Repair the broken PAP mapping reference.** Completed 29 September 2026. Removed the stale reference to the commented-out backup-status attrition table; retained the monitoring attrition, observability missingness, and Lee bounds references. Source: [PAP mapping](</home/ed/projects/overleaf/overleaf-takeup/tables/manual/pap-mapping.tex:94>).

- [x] **Apply the reviewed missingness wording and diagram.** Completed 29 September 2026. Added the approved main-text footnote explaining 252/126, revised the missingness-table sample description, and changed the PAP count to 252. Subsequently replaced the sample-construction paragraph with the detailed 1,204 individual / 1,203 paired / 252 missing / 1,141 usable accounting. Simplified the diagram to a branch from the 2,659-person endline sample to the 1,141-person observability sample, labeled “Individual observability measure available.” Retained the original closing table-note explanation as requested. PDF builds and diagram visually checked. The separate upstream 3,729 field-completion versus 3,678 clean-record discrepancy remains unresolved; this completion mark does not certify those totals. See the [missingness audit](observability-missingness-audit-2026-09-29/README.md).

- [x] **Finish the editorial pass.** Completed 29 September 2026 with the three approved edits in `experimental-design.tex`: lead the take-up paragraph with the distance finding, describe the 9,805-person sample as not enrolled in SMS, and describe the valuation exercise as a switching decision at a randomly assigned cash offer. Retained the broader willingness-to-accept wording in the introduction; other table-first openings were reviewed as optional polish rather than necessary corrections.

- [x] **Inspect the final PDF.** Completed 29 September 2026: visual overview of all 122 pages plus full-size checks of problem pages and the revised sample diagram. No undefined references. See the [inspection report](final-pdf-inspection-2026-09-29/README.md).

- [ ] **Resolve the remaining PDF layout findings.** Three oversized tables have notes colliding with page numbers: B15 (p. 64), L1 (p. 115), M1 (p. 118). Keep the appendix M heading with its table (p. 117), avoid the stranded PAP group heading (p. 96), and break the final overflowing formula (p. 122). Rebuild and inspect affected pages afterward.

- [x] **Apply the upstream endline count reconciliation.** Completed 29 September 2026 after approval. Updated the design, appendix, PAP, and diagram to 3,776 distinct reached respondents in the cleaned survey records, 98 who did not complete, and 3,678 completed respondents. The appendix identifies 1,019 SMS recipients excluded to obtain 2,659. Applied “Did not complete: 98” without the unsupported refusal/unable split. Synced the source diagram and frozen replication copy and refreshed its artifact-contract hash. The read-only audit reproduces every saved respondent/submission-key pair; no estimation changes required. See the [completion reconciliation](observability-missingness-audit-2026-09-29/README.md#follow-up-endline-completion-reconciliation).

## Already addressed — avoid duplicating work

- [x] Explain why structure is needed without a zero-observability arm.
- [x] Add a structural-section introduction.
- [x] Add descriptions of the previously missing alternative models. Review resolved by the user.
- [x] Discuss observability sensitivity in the structural section.
- [x] Update benchmark policy counts to 102/77 and the posterior draw count to 1,600.
- [x] Correct the distance-density caption's description of observed support.

These completion marks record the presence of the relevant manuscript changes, not independent certification of every underlying calculation.
