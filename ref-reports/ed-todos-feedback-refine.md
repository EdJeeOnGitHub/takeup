# Ed ToDos from Refine Feedback

Extracted from `ref-reports/edited-feedback-refine.md`. Excludes items whose `ToDo` is `AK`.

Paper source: the `.tex` paper files are in `../overleaf/overleaf-takeup`; the current draft is `ECM ReStud.tex`.

## Atomic commit map

All paper-text commits below are on the Overleaf branch `ed-refine-todos`. The commit diff is the coauthor-facing `.tex` diff for that ToDo.

- **ToDo 1, prior deworming**: Analysis `54fa50d` adds the generator script; `7705024` revises it to report the full Bracelet/Calendar/Ink/Control by Close/Far specification with prior-deworming covariates shown at the bottom; `7cfa8f3` removes generated RF table output from tracking so the table remains reproducible rather than source-controlled. Overleaf `b78969a` adds the Calendar-placebo/prior-deworming discussion, defines the baseline prior-deworming variables, and inputs `tables/prior-deworming-robustness`; `9baa834` refreshes the paper text/note/table for the full arm-by-distance specification.
- **ToDo 2, balance control means**: Analysis `f9d733e` changes balance table generation to report raw Control means from the relevant estimation sample. No `.tex` paper edit.
- **ToDo 4, \(\psi\) units**: Analysis `219dc73` labels generated WTP/structural output in USD. Overleaf `15fc692` changes the Section 4.2.2/Appendix WTP wording, WTP table note, and structural parameter table row so \(\psi\) is reported as USD.
- **ToDo 5 and 8, potential-distance/continuous-distance support**: Overleaf `9f5df01` clarifies that simulations condition on realized surveyed communities and rerun PoT selection, Close/Far assignment, targeted-community pairing, and feasibility screening. No separate analysis commit beyond the existing simulated-assignment inputs.
- **ToDo 6, \(\sigma_u\) identification**: Analysis `fe9c8d8` adds the prior/posterior plotting script and summary output; `7705024` updates the plot to use the paper's Canva palette. Overleaf `10cfa9c` adds the identification text and Appendix Figure `fig:sigma-u-prior-posterior`; `9baa834` refreshes the figure and normal-sum wording.
- **ToDo 10, Figure 3(b) curve ordering**: Overleaf `c013bdb` adds the clarifying sentence that the lower-noise curves remain above higher-noise curves at fixed cutoff. No analysis code change.
- **ToDo 11, HMC diagnostics**: Analysis `e7e2ac1` adds the reusable HMC diagnostics script. Overleaf `2a3ca83` replaces the loose main-table convergence sentence with checked max split \(\hat R\), ESS, and divergence stats; `9831b83` adds the corresponding diagnostics to the three appendix structural robustness table notes.
- **ToDo 12, E3 continuous-distance RI**: Analysis `41873f6` replaces free within-county distance shuffles with feasible counterfactual distance pools from the constrained assignment simulations. Overleaf `90a16d5` updates the E3 text/note and figure.
- **ToDo 17, Figure M1 parameters**: Overleaf `7df4781` updates Section M and the Figure M1 note to match the plotted Gaussian and bimodal distributions and clarify that \(V^\ast\) is the \(w^\ast\) cutoff scale.
- **ToDo 18, closest PoT and original-distance ITT check**: Overleaf `9195391` corrects the closest-PoT/candidate-center language. Analysis `fe26587` adds `--itt` reduced-form generation and merge checks for original distance assignment; `a22d7eb` fixes LaTeX escaping for small p-values in generated RF tables. Generated RF tables are intentionally untracked.
- **ToDo 24, terminology**: Overleaf `0dbc7bb` replaces "normal mixture"; `9baa834` updates the final wording to "normal-sum specification." No analysis code change.
- **ToDo 25, SMS control-mean SE and custom table style**: Analysis `b59c86b` changes the SMS control mean parenthetical to a clustered standard error and adds compact table postprocessing for SMS/heterogeneity table generation. Generated RF tables are intentionally untracked. No `.tex` paper edit.

---

## 1. Severe baseline imbalance in past deworming aligns with main effects

**ToDo**: Check if previous worming is predictive of takeup in the paper (my guess no based on double lasso), bound TEs if so.

Check definition of dewormed in the past year/dewormed in the past.

Within 1,995 check if predictive within that sample.

Use calendar as argument against:

>Before you move to your second step of bounding the treatment effect, you have a very powerful conceptual defense already in the paper: the Calendar arm. The baseline data shows that the Calendar arm has largely the same Far-Close imbalance in past deworming as the Bracelet and Ink arms. Yet, in the main results, the Calendar arm yields completely flat take-up across distance (it acts as a clean placebo). If the reviewer were correct that this specific baseline imbalance mechanically drives the take-up gap, the Calendar arm should have shown the same attenuated distance gradient as the public-signal arms. It did not.

**Implemented**: Checked the baseline prior-deworming variables. `treated_lgl` is prior/ever dewormed (`treated == "yes"`), and `treated_past_year` is treatment reported in the previous 1--12 months. The baseline survey IDs are not directly linkable to monitored take-up outcomes, so the predictive check uses cluster-level baseline shares merged to the full take-up sample. Added `scratch/prior-deworming-robustness-table.R`, which generates `presentations/rf-tables/main-specs/prior-deworming-robustness.tex`; copied the table to the paper tables directory and added Appendix Table `tab:prior-deworming-robustness` in `ECM ReStud.tex`. The generated table now reports the full Bracelet/Calendar/Ink/Control by Combined/Close/Far/Far--Close specification and shows the cluster shares ever dewormed and dewormed in the past year as covariates at the bottom. With those controls, the monitored Bracelet and Ink Far--Close interactions are 11.8 and 10.5 pp while Calendar is 2.7 pp; the endline self-reported specification shows the same qualitative pattern. Added text in `ECM ReStud.tex` defining the variables, using the Calendar placebo logic, and reporting this robustness check, so no bounding exercise is needed.

**Quote**:
> The baseline sample is broadly balanced across item/signal arms. Statistically significant differences are sparse relative to the number of reported comparisons and do not reveal a systematic pattern across treatment arms. In particular, the social image measures in Panel C are balanced across arms.

**Feedback**:
Table C1 appears to undercut the statement that baseline differences do not reveal a systematic pattern. The variable "Dewormed in the past year" shows a large and statistically significant imbalance: Control is higher than Bracelet/Ink in Close communities but lower in Far communities, so the implied Far-Close differences for Bracelet and Ink are much less negative than in Control. Because prior deworming is plausibly predictive of current take-up, this imbalance should be addressed when interpreting the distance-by-public-signal take-up interactions, even if it does not by itself invalidate the randomized design.

---

## 2. Mathematically impossible sample means in balance tables

**When**: 27-06-2026

**ToDo**: Rerun balance with control mean fixed for held out county strata.

**Implemented**: Updated `R/balance/functions.R` so displayed Control means in treatment/distance balance comparisons are raw means from the relevant estimation sample, not the omitted-county fixed-effect intercept. Updated `balance.R` to pass the underlying data into `create_balance_comparisons()` for the main, appendix, attrition, monitoring-attrition, and SMS-enrollment balance paths. Re-ran `Rscript balance.R --main --output-path=temp-data` and `Rscript balance.R --orig --output-path=temp-data` in the devcontainer; verified `temp-data/comp_takeup_balance_tidy_df.csv` and `temp-data/comp_balance_tidy_df.csv` now report raw Control means for the cited variables, e.g. Age Control-Close \(40.3\), Control-Far \(40.2\), Female Control-Close \(0.565\), Control-Far \(0.562\), Phone owner Control-Close \(0.765\), Control-Far \(0.728\).

**Quote**:
> | Age | 39.8 | 38.126 (0.758) | -0.714 [0.427] | -0.919 [0.286] | 0.232 [0.813] | 1.151 [0.208] | 0.5221 | 37.995 (0.838) | -0.804 [0.403] | -0.861 [0.324] | -0.798 [0.419] | 0.063 [0.945] | 0.7592 | 0.7776 |

**Feedback**:
In Table B1, the 'Sample mean' values for several variables fall completely outside the range implied by the 'Control mean' and treatment difference columns. For instance, the Sample mean for Age is 39.8, but combining the Control means with the corresponding treatment differences suggests all subgroup means correctly aggregate to a range between 37.1 and 38.4. The same pattern appears for the 'Female' and 'Phone owner' variables, and casades across Tables B2 and B12. This discrepancy apparent mechanically arises because the 'Control mean' column reports the intercept from the balance regression, representing the mean for the omitted county stratum rather than the raw overall Control mean. Presenting this omitted-stratum intercept as the generic 'Control mean' adjacent to the overall raw 'Sample mean' creates an apparent contradiction that is confusing for readers trying to deduce baseline characteristics from the table. To resolve this, consider either reporting raw Control means, employing fixed-effect adjusted means centered at the sample average of the fixed effects, or explicitly noting in the table notes that the 'Control mean' corresponds to the regression intercept for the omitted stratum.

---

## 3. Post-treatment selection in gift-choice test

**ToDo**: AK - not for us. Ignore.

**Quote**:
> One might also worry that private value depends on local item prevalence. We therefore examine gift choices among non-dewormers in Bracelet and Calendar communities who did not already receive an item through deworming. Table B24 shows item prevalence effects: when the calendar is locally common, non-dewormers are more likely to choose the bracelet than in Control communities ( $+17.3 \mathrm{pp}, p<0.01$ ). When the bracelet is locally common, non-dewormers are less likely to choose the bracelet than in Control communities $(-9.7 \mathrm{pp}, p<0.05)$.
>
> If marginal adults in Far Bracelet communities had a higher private valuation for the bracelet, bracelet choice should increase in Far Bracelet communities relative to comparable adults in Far Control communities. It does not: relative to Control-Far, bracelet choice in Bracelet-Far is lower by 8.9 percentage points ( $p<0.10$; Table B24). The Far-Close difference in the Bracelet coefficient is also small and statistically insignificant (+1.6 pp, s.e. 6.3 pp ).

**Feedback**:
In Section 4.2.2, the item-prevalence and Far Bracelet comparisons condition on being a non-dewormer, which is itself affected by the item/signal assignment. Since higher-bracelet-value adults in Far Bracelet communities may be precisely the people induced to deworm and then excluded from this sample, the -8.9 pp estimate does not by itself rule out high bracelet valuation among the relevant marginal compliers.

---

## 4. Unit error in relative valuation scale and reported estimates

**ToDo**: Change B23 to report \psi in USD.

**Implemented**: Updated `ECM ReStud.tex` notes and Appendix WTP likelihood text to define \(\psi\) in USD; updated `quick_param_postprocess.R` and `presentations/tables.Rmd` so generated structural/WTP tables label the valuation difference in USD.

**Quote**:
> Let
>
> $$
> g_{i}^{\mathrm{val}} \equiv U_{c, i}-U_{b, i}
> $$
>
> denote respondent $i$ 's latent private value of the calendar relative to the bracelet, measured in the monetary units used in the elicitation. We assume
>
> $$
> g_{i}^{\mathrm{val}} \sim N(\psi, 1),
> $$
>
> where $\psi$ is the mean calendar minus bracelet value and the variance is normalized for scale.

**Feedback**:
There appears to be a unit inconsistency in the relative valuation model/reporting. Appendix J says $g_i^{\mathrm{val}}$ is measured in the monetary units used in the elicitation and fixes $g_i^{\mathrm{val}} \sim N(\psi,1)$, while Table B23 reports $\psi=0.472$ as KSh and also reports nondegenerate choice probabilities for 50 KSh offers. Furthermore, the text in Section 4.2.2 explicitly interprets the estimated mean gap as US$0.47, which implies a value about 100 times larger than 0.472 KSh. If $m_i$ entered the likelihood as raw KSh, a 50 KSh offer would be a 50-standard-deviation shift under the stated normalization; the reported probabilities instead suggest that the offer amounts were scaled, such as in dollars or hundreds of KSh. It is not fully clear whether $\psi$ is a monetary valuation or a normalized probit-index parameter, which matters for interpreting the valuation estimates and the mapping through $\rho$. The monetary scale used for $m_i$ and $\psi$ needs to be clarified.

---

## 5. Appendix E potential-distance check is underspecified

**ToDo**: Flesh this out and explain, make clear how we generate.

**Implemented**: Expanded Section 3.1 and Appendix E text in `ECM ReStud.tex` to state that simulations hold surveyed communities fixed, rerun PoT selection/Close-Far assignment/targeted-community pairing/feasibility screening, assign potential distances under simulated feasible pairings, and compare equally weighted community-draw distributions.

**Quote**:
> A potential concern is that the acceptance-rejection algorithm may have introduced systematic differences in the set of feasible distances for communities assigned to Close versus Far, complicating the interpretation of realized continuous distance. Following Borusyak and Hull (2023), we assess this by simulating the realized selection and assignment mechanism. Holding the set of surveyed communities fixed, we re-run the full PoT selection and assignment algorithm 500 times and compute the centroid-to-PoT distances that would be induced under these counterfactual PoT-selection and assignment re-draws. We then compare the resulting potential distance distributions for communities assigned to Close and Far in the realized experiment. Figure E2 shows that these distributions are statistically indistinguishable, providing no evidence that the constrained algorithm generated systematic differences in potential distance profiles across realized assignment groups.

**Feedback**:
Appendix E’s potential-distance diagnostic is difficult to reconstruct from the description. In particular, it does not specify how each fixed surveyed community is assigned a counterfactual PoT distance when a simulated re-run selects a different set of PoT centers or targeted-community pairings, nor does it state the unit or weighting used for the Kolmogorov-Smirnov comparison.

---

## 6. Missing identification source for noise parameter

**ToDo**: it's the first two (curvature in takeup across cells + nonlinear normal specification). Plot prior and posterior for \sigma_u.

**Implemented**: Updated the identification paragraph in `ECM ReStud.tex` to state that \(\sigma_u\) is disciplined by curvature in arm-by-distance take-up cells plus the nonlinear normal-sum signal-extraction restriction. Added `scratch/sigma-u-prior-posterior-plot.R`, which reads only the `u_sd` Stan variable, generates `presentations/figures/sigma-u-prior-posterior.pdf`, and writes `presentations/figures/sigma-u-prior-posterior-summary.csv`. Copied the figure to the paper figures directory and added Appendix Figure `fig:sigma-u-prior-posterior` in the priors subsection. The generated summary gives prior median \(0.374\) and posterior median \(0.278\).

**Quote**:
> Because latent utility is normalized by $\operatorname{Var}\left(v_{i}\right)=1, \sigma_{u}$ is interpreted as idiosyncratic decision noise relative to intrinsic motivation. The joint likelihood pins down the implied social image return, $\lambda p_{\text {Observed }}(z, d) \Delta\left(w^{*}\right)$, and how it changes with distance. These are the objects needed for the unitless social multiplier.
>
> Remaining parameters. Conditional on the relative valuation parameter $\psi$, the noise parameter $\sigma_{u}$, and the observability parameters $\beta_{z}^{O}$ and $\gamma_{z}^{O}$, the remaining take-up parameters are identified through the take-up likelihood.

**Feedback**:
There seems to be an issue with the identification account for $\sigma_u$: the subsection explains how $\sigma_u$ affects the type-inference term and the multiplier, but it does not identify the empirical variation or parametric restriction that separately disciplines it. Because the later priors discussion notes that $\sigma_u$ requires regularization, the current identification paragraph leaves unclear whether this parameter is primarily identified by the nonlinear normal specification, by curvature in take-up across distance/observability cells, or materially by the prior.

---

## 8. Section 3.1 continuous-distance exogeneity needs support

**ToDo**: Take Refine's sentence and add to 3.1

**Implemented**: Added the conditioning-set language to Section 3.1 in `ECM ReStud.tex`, clarifying that continuous realized distance is interpreted conditional on the fixed surveyed-community frame and constrained assignment protocol.

**Quote**:
> Distance as a continuous measure. Our primary reduced form distance treatment is the Close/Far assignment. We use realized community-centroid-to-PoT distance only in reduced form robustness checks, where it provides a continuous measure of the travel distance variation induced by the design. In the structural model, centroid-to-PoT distance is the main distance input, since it is the randomized spatial margin and a salient component of access costs. Specifically, we calculate the distance from each community's centroid to its assigned finalized PoT. All realized distance measures use the finalized PoT location at which treatment was delivered.
>
> Because PoTs are selected using a geographically constrained procedure, the induced assignment mechanism could in principle depend on local geography, such as school density. Following Borusyak and Hull (2023), we re-run the full selection and assignment algorithm 500 times to assess whether communities assigned to Close and Far in the realized experiment had similar potential distance distributions under counterfactual re-draws. These simulations provide no evidence of systematic differences in potential centroid-to-PoT distances across the realized Close and Far groups (Figure E2). We further verify that realized centroid-to-PoT distance is uncorrelated with baseline covariates (Figure E3). Taken together, these checks support interpreting realized community-centroid-to-PoT distance as a continuous measure of the as-good-as-random travel distance variation induced by the design.

**Feedback**:
The discussion of realized centroid-to-PoT distance could more sharply distinguish the randomized Close/Far assignment from the exact continuous distance realized after community selection, finalized PoT locations, and re-randomizations. The simulations and covariate checks support the interpretation, but it remains somewhat unclear what conditioning set justifies treating exact realized distance as as-good-as-random in the structural model.

---

## 10. Impossible crossing of curves in Figure 3(b)

**ToDo**: This seems wrong, the curves don't cross?

**Implemented**: Checked the current Figure 3(b) data in `temp-data/delta-vstar-endog-mu.csv`. Interpolating the highlighted \(\sigma_u=0.5,1,1.5\) curves on their common \(w^\ast\) range shows no crossing: the lower-noise curve remains above the higher-noise curve. Added a clarifying sentence to `ECM ReStud.tex` stating that Panel (b) varies \(\mu(d)\) with distance and does not reverse the ordering of \(\Delta(w^\ast)\) across \(\sigma_u\).

**Quote**:
> Panel (b) allows observability to fall with distance, $\mu^{\prime}(d)<0$. If the decline in $\mu(d)$ dominates the rise in $\Delta$, then the social image return can fall with distance: at large $d$, only highly motivated individuals deworm, but their action is less visible, so the reputational return $\mu(d) \Delta\left(w^{*}(d)\right)$ is small.

**Feedback**:
Figure 3(b) appears to show the social-image-return curves crossing so that, for some positive $w^*$, a smaller idiosyncratic-noise dispersion such as $\sigma_u=0.5$ yields a lower $\mu(d)\Delta(w^*)$ than a larger dispersion such as $\sigma_u=1.5$. Under the stated normal model and a common declining observability schedule, $\Delta(w^*)$ should fall as $\sigma_u$ rises at any fixed cutoff, and the implicit equilibrium mapping should preserve this ordering when $\mu'(d)\le 0$. The lower-$\sigma_u$ curve should therefore remain weakly above the higher-$\sigma_u$ curve; the apparent crossing suggests a possible inconsistency in how Panel (b) maps cutoff, distance, and observability.

---

## 11. Outdated $\hat{R}$ convergence threshold for HMC

**ToDo**: Check whatever latest value is and verify we reach.

**Implemented**: Replaced the loose \(\hat R<1.1\) statement in `ECM ReStud.tex` with checked diagnostics for the reported parameters in `tab:struct-model-params`: max split \(\hat R=1.011\), min bulk ESS 465, min tail ESS 602, and zero divergent transitions. Added reusable script `scratch/stan-hmc-diagnostics.R` for model-specific diagnostics with Stan variable filtering. Running it on the no-outlier robustness fit (`fit105`, `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS`) gives max split \(\hat R=1.018\), min bulk ESS 355, min tail ESS 309, and zero divergent transitions. Main-fit downstream tables inherit the main structural diagnostics rather than needing repeated table-note Rhats. The individual-distance appendix model diagnostics were computed on Midway: community-FP has max split \(\hat R=1.022\), min bulk ESS 221, min tail ESS 212, and zero divergences; full-info individual-FP has max split \(\hat R=1.022\), min bulk ESS 142, min tail ESS 132, and two divergences. Structural-table audit:

| Paper table / figure | Table input in `.tex` | Producing model or fit | Found locally? | \(\hat R\) stats added? |
|---|---|---|---|---|
| `tab:wtp-summ-table` | `tables/wtp-summ-table` | Joint baseline structural fit, WTP block | Yes: table and baseline fit outputs present | Covered by main structural diagnostics |
| `tab:struct-model-params` | `tables/all-structural-parameters-table-fit-1-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB` | `dist_fit95_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB`, chains 1--4 | Yes: raw fit outputs present | Yes; note reports max \(\hat R\), min ESS, and divergences for reported parameters |
| `tab:beliefs-table` | `tables/fit105/fob-beliefs-table` | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP`, fit 105 postprocessed belief ATEs | Yes: table, postprocessed RDS, and raw fit outputs present | Covered by main structural diagnostics |
| `tab:struct_overall_effects` | `tables/fit105/struct-overall-te-table` | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP`, fit 105 postprocessed ATEs | Yes: table, postprocessed RDS, and raw fit outputs present | Covered by main structural diagnostics |
| `tab:signaling_private` | `tables/fit105/private-signal-te-table` | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP`, fit 105 postprocessed private/social-image ATEs | Yes: table, postprocessed RDS, and raw fit outputs present | Covered by main structural diagnostics |
| `tab:optim-post` | `tables/fit105/optim-summ-table` | Policymaker optimization using posterior draws from `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP`, fit 105 | Yes: table and optimization summaries present; raw fit outputs present | Covered by main structural diagnostics |
| `fig:sigma-u-prior-posterior` | `figures/sigma-u-prior-posterior.pdf` | Prior and posterior draws of `u_sd` from `dist_prior95/dist_fit95_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB` | Yes: raw prior/posterior outputs present | Not a convergence-diagnostics figure |
| `tab:indiv-community-struct` | `tables/fit105/indiv-dist-community-fp-indiv-vis-robust-struct-overall-te-table` | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS`, fit 105 postprocessed robustness ATEs | Table and postprocessed RDS present locally; raw model-suffixed fit outputs found on Midway | Computed separately: max \(\hat R=1.022\), min bulk ESS 221, min tail ESS 212, zero divergences |
| `tab:indiv-indiv-struct` | `tables/fit105/indiv-dist-indiv-fp-robust-struct-overall-te-table` | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP`, fit 105 postprocessed robustness ATEs | Table and postprocessed RDS present locally; raw model-suffixed fit outputs found on Midway | Computed separately: max \(\hat R=1.022\), min bulk ESS 142, min tail ESS 132, two divergences |
| `tab:nooutlier-struct` | `tables/fit105/struct-robustness-nooutliers-overall-te-table.tex` | `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS`, fit 105 postprocessed robustness ATEs | Yes: table, postprocessed RDS, and raw no-outlier fit outputs present | Computed separately: max \(\hat R=1.018\), min bulk ESS 355, min tail ESS 309, zero divergences |
| `tab:sigmau-robustness` | inline table in `.tex` | Posterior summaries of \(\sigma_u^2\) for main and full-information specifications | Main raw fit outputs present locally; full-information model outputs found on Midway | Main side covered by main diagnostics; full-information side covered by `tab:indiv-indiv-struct` diagnostics |
| `tab:optim-post-robust` | `tables/fit105/optim-summ-robust-table` | Policymaker optimization using posterior draws from `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP`, fit 105 | Yes: table and optimization summaries present; raw fit outputs present | Covered by main structural diagnostics |

**Quote**:
> The model is estimated in Stan using Hamiltonian Monte Carlo with four chains, 400 warm up draws, and 400 saved draws per chain. Sampler diagnostics show no divergent transitions and split $\hat{R}<1.1$ for all parameters.

**Feedback**:
The HMC diagnostics are reported using a relatively loose convergence threshold. For a Bayesian structural model that feeds into decomposition and policy counterfactuals, $\hat{R}<1.1$ alone does not fully rule out mixing problems, especially with 400 saved draws per chain. If the actual values are all comfortably below stricter thresholds and effective sample sizes are adequate, this is only a reporting issue; otherwise, the posterior summaries may require longer runs or additional convergence evidence.

---

## 12. Appendix E.1 permutation test may not match the design

**ToDo**: Re-run figure E.3 using the algorithm.

**Implemented**: Replaced the free distance-shuffle RI in `balance.R` with placebo continuous-distance draws from `data/simulated-counterfactual-treatment-assignment.csv`, sampling feasible counterfactual PoT-community distances for each observed community from the constrained assignment simulations. Re-ran `Rscript balance.R --orig --fit-ri --output-path=temp-data` in the devcontainer, regenerated `temp-data/balance-cts-dist-ri-fe.rds`, `temp-data/balance-ri-fe.rds`, and `temp-data/dist-ri-plot.pdf`, and copied the updated E3 figure to `../overleaf/overleaf-takeup/figures/dist-ri-plot.pdf`. The minimum displayed covariate p-value is \(0.124\). Updated Appendix E.1 in `ECM ReStud.tex` to describe the simulation-based placebo draws rather than within-county permutations.

**Quote**:
> As an additional check that our continuous distance measure is orthogonal to pretreatment characteristics, we regress each baseline covariate on the distance from the community centroid to its assigned point of treatment (PoT), controlling for county fixed effects and clustering standard errors at the community level. We then permute community-level distances across communities within county strata 500 times and recompute the associated test statistic. Figure E3 plots the resulting randomization distribution for each covariate under the null of no association. The vertical line marks the realized value, and permutation $p$-values based on the 500 draws are reported in the upper right-hand corner of each panel. For all covariates shown, the $p$-value exceeds 0.10 , indicating that realized centroid-to-PoT distance is not systematically related to observable characteristics.

**Feedback**:
Appendix E.1’s permutation exercise is useful as a descriptive balance check, but it is not clear that freely permuting realized community-level distances within county reproduces the constrained PoT-selection and community-assignment mechanism. If the reported values are intended as design-based randomization-inference $p$-values, the exchangeability assumption behind this permutation scheme needs to be justified.

---

## 17. Contradiction between text and Figure M1

**ToDo**: Add correct numerical figures to text.

**Implemented**: Updated Section M in `ECM ReStud.tex` to describe the simulation parameters shown in Figure M1: the plotted Gaussian benchmark is approximately \(N(0.4,0.2^2)\), and the bimodal case is approximately \(0.25N(-2,0.2^2)+0.75N(0.4,0.2^2)\). The figure note now also clarifies that the plotted \(V^\ast\) label corresponds to the cutoff/type index scale \(w^\ast\) used in the text.

**Quote**:
> ![](/documents/b1b14a42-0216-4727-8e8b-e1392684c111/images/image_021.jpg)
> Figure M1: Simulation under Gaussian and Bimodal Type Distributions
>
> Notes: Panel (a) plots the assumed densities for $w$: Gaussian in red and bimodal in blue.

**Feedback**:
There appears to be a disconnect between the simulation parameters described in the text of Section M and those plotted in Figure M1. The text specifies a standard normal $w \sim N(0,1)$ and an "equally weighted bimodal mixture of $N(-\eta, 1)$ and $N(\eta, 1)$." However, Panel (a) of Figure M1 displays a highly asymmetric bimodal distribution (modes at approximately 0.4 and -2.0, with differing peak heights) and a Gaussian benchmark with a peak density of 2.0, which implies a standard deviation of approximately 0.2 rather than 1. Additionally, the axes in panels (a) and (c) are labeled $v^*$ and $V^*$, rather than $w$ or $w^*$ as denoted in the text and caption. Updating the text to accurately reflect the simulation parameters plotted in the figure would resolve this confusion.

---

## 18. Apparent contradiction regarding closest PoT

**ToDo**: Ed check whether each community is closest to it's own PoT.

**Implemented**: Clarified in Appendix E of `ECM ReStud.tex` that the 0.5 km separation guarantee applies to candidate PoT centers used by the spatial algorithm, while finalized implemented PoT locations are a separate object. Verified the realized implementation distances in `data/takeup_village_pot_dist.RData`: `vill.pot.dist` has all 144 selected communities with `dist.to.own.pot <= closest.other`, so every assigned finalized PoT is closest. The “two” issue is not closestness: in `verify.vill.pot.dist`, two far Kakamega clusters (1152 and 1236) have positive nearest-alternative margins slightly below 0.5 km (about 477 m and 471 m). Updated the paper to report this correctly.

**Additional check**: Added an original-distance-assignment ITT robustness to `scripts/reduced-form/bootstrap.R --itt`, using the existing reduced-form table writer. The code replaces the current Close/Far distance cell with `assigned.dist.cat` from `data/rct_targetable_schools_2.0.rds`, writes merge diagnostics to `temp-data/original-distance-merge-checks.csv`, switch counts to `temp-data/original-distance-switch-counts.csv`, tidy output to `temp-data/tidy-rf-tes/original-distance-itt-tidy-tes.csv`, and tables to `presentations/rf-tables/main-specs/rf_original_distance_itt_tbl.tex` and `_weird_order.tex`. The merge checks verify 9,805 rows, 144 clusters, zero missing original assignments, zero missing expected-distance controls, zero treatment mismatches against monitored and processed randomization data, zero current-distance-group mismatches, and PoT distance differences only at floating-point tolerance. The original-vs-current distance switch count is 64 close/close, 10 close/far, 16 far/close, and 54 far/far. The original-assignment far--close public-signal estimate is 0.124 (SE 0.036, \(p=0.001\)). Also added original-distance ITT observability checks for first-order and second-order beliefs, writing `temp-data/tidy-rf-tes/original-distance-itt-fob-tidy-tes.csv`, `temp-data/tidy-rf-tes/original-distance-itt-sob-tidy-tes.csv`, `presentations/rf-tables/main-specs/rf_original_distance_itt_fob_tbl.tex`, and `presentations/rf-tables/main-specs/rf_original_distance_itt_sob_tbl.tex`. Observability merge checks in `temp-data/original-distance-observability-merge-checks.csv` verify 1,141 rows, 144 clusters, zero missing original assignments, zero missing expected-distance controls, zero treatment/current-distance mismatches against endline source data, and identical PoT distance fields. Relative to current-distance specifications, the first-order-belief signal far--close estimate is stable at 0.132 to 0.136, while bracelet--calendar increases from 0.148 to 0.242; for second-order beliefs, signal increases from 0.085 to 0.120 and bracelet--calendar from 0.123 to 0.287. Added a standalone side-by-side comparison wrapper, `presentations/original-distance-itt-comparison-note.tex`, and compiled `presentations/original-distance-itt-comparison-note.pdf`, showing current-distance and original-distance-ITT tables for take-up and first-order observability.

**Quote**:
> Study communities are selected from the remaining, non-overlapping portion of the 2.5 km catchment. As a result, any selected study community is at least 0.5 km closer to its assigned PoT than to the nearest alternative PoT.

**Feedback**:
The 0.5 km separation statement appears to refer to the candidate PoT centers used in the spatial selection algorithm, but the later statement that in two clusters the assigned finalized PoT was not the closest deworming site creates a tension as written. The distinction between candidate PoT centers and finalized PoT locations should be made explicit at the point where the 0.5 km guarantee is stated.

---

## 24. Incorrect statistical terminology: normal mixture

**ToDo**: Ed edit the term.

**Implemented**: Replaced “normal mixture” with “normal-sum specification” in `ECM ReStud.tex`.

**Quote**:
> We normalize $v_{i} \sim N(0,1)$ and assume $u_{i} \sim N\left(0, \sigma_{u}^{2}\right)$, so $w_{i}=v_{i}+u_{i} \sim N\left(0,1+\sigma_{u}^{2}\right)$. For the normal mixture, the type inference term is

**Feedback**:
The phrase “normal mixture” is statistically imprecise here. Since $w_i=v_i+u_i$ is defined as the sum of independent normal variables and is immediately stated to be normal, this is a normal-sum/convolution specification rather than a mixture distribution.

---

## 25. Possible error in standard error reporting in Table B21

**ToDo**: Swap to standard error.

**Implemented**: Updated `scripts/reduced-form/bootstrap.R` so the SMS control mean parenthetical is a clustered standard error from an intercept-only regression, not the binary outcome SD. Regenerated `presentations/rf-tables/main-specs/incentive-control-sms-te.tex`; the control mean now reports `(0.032)` instead of `(0.470)`. Also added compact table postprocessing so the SMS and heterogeneity table inputs match the manually edited ECM ReStud layout/style, leaving notes/significance text in the paper. Fixed the `--heterogeneity` path to construct the cleaned WTP `second_choice`/county fields before WTP checks, so `Rscript scripts/reduced-form/bootstrap.R --heterogeneity --table-output-path=presentations/rf-tables/main-specs` now exits cleanly after regenerating `het-tes.tex`.

**Quote**:
> | Social Info | $0.131^{* * *}$ |
> |  | $(0.028)$ |
> | Control mean | 0.330 |
> |  | $(0.470)$ |
> | Observations | 1,705 |
>
> Notes: This table reports estimates from a linear probability model (OLS) with a binary outcome for deworming take-up. The regressors are indicators for receiving a Social Information SMS or a Reminder Only SMS. The omitted group is phone owners who received no SMS. Because the Reminder treatment was implemented only in Control communities, the analysis sample is restricted to phone owners in the Control arm, not the full SMS treatment sample ( $N=3,022$ ). All regressions include county fixed effects and the LASSO selected controls for sex and age, and expected distance as a design control. Bootstrapped standard errors, clustered at the community level, are in parentheses. ${ }^{*} p<0.10,{ }^{* *} p<0.05,{ }^{* * *} p<0.01$.

**Feedback**:
In Table B21, the parenthetical value below the Control mean appears to be the control-group standard deviation of the binary take-up outcome, not a bootstrapped clustered standard error. This conflicts with the note stating that values in parentheses are standard errors.
