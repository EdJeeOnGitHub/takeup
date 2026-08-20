# ToDo .tex Review Checklist

Use this as a checklist for reviewing every ToDo-related `.tex` change that landed in the paper or in paper table inputs.

Paper repo: `../overleaf/overleaf-takeup`  
Landed Overleaf branch: `overleaf/master`  
Main paper file: `ECM ReStud.tex`

To inspect an exact commit diff from the Overleaf repo:

```sh
git -C ../overleaf/overleaf-takeup show <commit> -- 'ECM ReStud.tex'
git -C ../overleaf/overleaf-takeup show <commit> -- <table-or-rf-table-path>
```

## Review Items

### ToDo 1: Prior Deworming Robustness

- [ ] Review `ECM ReStud.tex`.
  - Landed commit: `0d33595 todo1: add monitored prior deworming robustness`
  - Main text target: the Calendar/placebo paragraph around Appendix Table `\ref{tab:prior-deworming-robustness}`.
  - Check that the text only discusses monitored take-up and cluster-level baseline prior-deworming controls.
  - Check that it does not discuss a self-reported deworming outcome.

- [ ] Review `tables/prior-deworming-robustness.tex`.
  - Landed commit: `0d33595 todo1: add monitored prior deworming robustness`
  - Check that the table is the full Bracelet/Calendar/Ink/Control by Close/Far specification.
  - Check that prior-deworming covariates appear at the bottom of the table.

### ToDo 2: Balance Table Control Means

- [ ] Review `tables/stacked-sample-balance-table.tex`.
  - Landed commit: `4640975 todo2: update balance table control means`
  - Check that Control means are raw Control means for the relevant sample, not held-out-strata means.

- [ ] Review `tables/stacked-sample-dist-balance-table.tex`.
  - Landed commit: `4640975 todo2: update balance table control means`
  - Check that the distance-split balance table uses the corrected Control mean logic.

- [ ] No `ECM ReStud.tex` prose change was needed for this ToDo.

### ToDo 3: AK ToDo

- [ ] No paper `.tex` change for us to review.

### ToDo 4: \(\psi\) Units

- [ ] Review `ECM ReStud.tex`.
  - Landed commit: `7fe05b8 todo4: clarify psi valuation units`
  - Check the Section 4.2.2 and appendix WTP/structural wording for `\psi`.
  - Confirm the text says the valuation difference is in USD/monetary units from the elicitation likelihood.

- [ ] Review `tables/wtp-summ-table.tex`.
  - Landed commit: `0de254e todo4: label WTP valuation in USD`
  - Check the row label `Valuation difference (USD), mean`.

### ToDos 5 and 8: Potential-Distance and Continuous-Distance Support

- [ ] Review `ECM ReStud.tex`.
  - Landed text commit: `905c279 todo5: clarify potential distance simulations`
  - Landed figure-refresh commit: `ac9a5ac todo5: refresh counterfactual distance plot`
  - Check the site-selection appendix language describing the counterfactual assignment simulation.
  - Check that the figure block points to `misc-figures/counterfactual-assignment-density-close-far-dens.pdf`.
  - Check that the text says the simulation conditions on realized surveyed communities and reruns feasible PoT/community assignment steps.

### ToDo 6: Missing Identification / \(\sigma_u\)

- [ ] Review `ECM ReStud.tex`.
  - Landed commit: `68712c0 todo6: add sigma_u identification figure`
  - Later wording cleanup: `30215bd todo18 todo24: clean up PoT and normal-sum wording`
  - Check the structural-model identification wording for `\sigma_u`.
  - Check the priors appendix wording around Figure `\ref{fig:sigma-u-prior-posterior}`.
  - Check that the wording discusses the prior and posterior comparison without sounding defensive.

- [ ] Review the figure inclusion in `ECM ReStud.tex`.
  - Figure path: `figures/sigma-u-prior-posterior.pdf`
  - Label: `fig:sigma-u-prior-posterior`

### ToDo 10: Figure 3(b) Curve Ordering

- [ ] No final reader-facing `.tex` change landed on `overleaf/master`.
  - We reviewed this and decided the reviewer concern was not a real issue.
  - Branch-only review commit, not part of final master: `c013bdb todo10: clarify delta curve ordering`.

### ToDo 11: HMC Diagnostics

- [ ] Review the main structural table note in `ECM ReStud.tex`.
  - Landed commit: `9126723 todo11: report structural HMC diagnostics`
  - Check that the note reports divergent transitions, max split `\hat R`, minimum bulk ESS, and minimum tail ESS.
  - Check that it does not discuss Stan CSV files or implementation details that should not appear in a journal submission.

- [ ] Review the appendix structural robustness table notes in `ECM ReStud.tex`.
  - Landed commit: `4d1722e todo11: add appendix HMC diagnostics`
  - Check the notes for the household-distance private-cost table.
  - Check the notes for the household-distance private-cost and social-image table.
  - Check the notes for the no-outlier structural robustness table.

### ToDo 12: Appendix E3 Continuous-Distance RI

- [ ] Review `ECM ReStud.tex`.
  - Landed commit: `da36a88 todo12: update continuous distance RI note`
  - Check the Appendix E3 text and figure note.
  - Confirm it describes feasible counterfactual distance pools from the constrained assignment simulations.

- [ ] Review the figure inclusion in `ECM ReStud.tex`.
  - Figure path: `figures/dist-ri-plot.pdf`

### ToDo 17: Figure M1 Parameters

- [ ] Review `ECM ReStud.tex`.
  - Landed commit: `b05d100 todo17: clarify figure M1 parameters`
  - Check the Section M text and Figure M1 note.
  - Confirm the text matches the plotted Gaussian and bimodal distributions.
  - Confirm `type index` wording has no hyphen.

### ToDo 18: Closest PoT and Original-Distance ITT

- [ ] Review the closest-PoT wording in `ECM ReStud.tex`.
  - Landed commits:
    - `30215bd todo18 todo24: clean up PoT and normal-sum wording`
    - `206d986 todo18: align main PoT closestness wording`
  - Check the main text and site-selection appendix.
  - Confirm the wording says each selected community's assigned finalized PoT is the closest finalized deworming site.

- [ ] Review the original-distance ITT table insertions in `ECM ReStud.tex`.
  - Landed commit: `d4e36e4 todo18: add original distance ITT reduced form tables`
  - Check that these appear in the robustness section for reduced-form models.
  - Check that the paper includes take-up and FOB observability only.
  - Confirm SOB observability is not included.

- [ ] Review `rf-tables/main-specs/rf_original_distance_itt_tbl.tex`.
  - Landed commit: `d4e36e4 todo18: add original distance ITT reduced form tables`
  - This is the take-up ITT table using the original distance assignment.

- [ ] Review `rf-tables/main-specs/rf_original_distance_itt_fob_tbl.tex`.
  - Landed commit: `d4e36e4 todo18: add original distance ITT reduced form tables`
  - This is the FOB observability ITT table using the original distance assignment.

### ToDo 24: Normal-Sum Wording

- [ ] Review `ECM ReStud.tex`.
  - Landed commit: `30215bd todo18 todo24: clean up PoT and normal-sum wording`
  - Check the structural-model wording around `w_i=v_i+u_i`.
  - Confirm the final wording uses `normal-sum specification`.
  - Confirm it does not use `additive normal specification`.

### ToDo 25: SMS Control Mean Standard Error

- [ ] Review `rf-tables/main-specs/incentive-control-sms-te.tex`.
  - Landed commit: `ade3758 todo25: update SMS control mean standard error`
  - Check that the Control row parenthetical is the clustered standard error, not the Control standard deviation.

- [ ] No `ECM ReStud.tex` prose change was needed for this ToDo.

## Review-Packet .tex Artifacts

These are not paper edits, but they are `.tex` files created for our review workflow.

- [ ] Review `presentations/todo-coauthor-review.tex`.
  - Purpose: coauthor-facing summary packet.
  - Current title should be `ToDo Reviews`.
  - It should not include the sentence `This packet is meant for coauthor review`.
  - It should not include the removed self-reported deworming sentence.

- [ ] Review `presentations/todo-change-review.tex`.
  - Purpose: more technical change-review packet with tables, figures, and diffs.

## Items With No Final .tex Change

- [ ] ToDo 3: AK-owned item, no final `.tex` change for us.
- [ ] ToDo 10: intentionally ignored after review, no final `.tex` change on `overleaf/master`.
- [ ] Sigma-u plot implementation details: no journal-facing `.tex` should mention Stan CSV loading or local compute details.
