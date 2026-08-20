# The Social Multiplier of Access Costs: Experimental Evidence from Community Deworming in Kenya*

**Date**: 5/26/2026, 10:26:20 PM
**Domain**: social_sciences/economics
**Taxonomy**: academic/research_paper
**Filter**: Active comments

---

## Overall Feedback

Here are some observations following a careful read of the manuscript.

**Experimental anchoring of the direct travel-cost effect**

The headline contribution features an estimated social-image multiplier of 0.8–1.6. Defining this multiplier relies on isolating the ratio of the total distance effect to the direct distance effect that would obtain if the social image channel were completely shut down (Sections 5.2 and 6.1.2).

A central question for the structural decomposition is how the direct travel-cost denominator is experimentally anchored. In the experiment, Control communities do not represent a pure private-information environment; Table 1 shows that reported observability remains high in Close communities and falls with distance. Consequently, the direct travel-cost slope appears to be recovered via the joint structural restrictions on distance costs, observability, the social-image weight, the type distribution, and idiosyncratic noise, rather than from a transparent experimental setup with zero observability. Presenting explicit evidence on which moments independently pin down the direct distance slope versus the reputational feedback will substantiate the multiplier claim.


**Mapping the observability measure to social image payoffs**

Section 6 maps $p_{Observed}(z,d)$ to the probability that an individual's participation is known in the community, allowing the observability schedule to directly identify the social-image component. The primary empirical input is the respondent's share of recognized peers receiving a definite Yes/No answer.

An apparent empirical tension is that the structural model treats this measure as the exact object required for social-image payoffs, yet the measurement captures a noisier construct. Table B15 shows that public signals increase definite classifications while simultaneously lowering conditional accuracy in certain cells. Furthermore, the structural posterior does not presently propagate either misclassification or the differential missingness in the observability module documented in Tables B4–B6. Showcasing that the multiplier estimates survive alternative observability inputs—such as unconditional correct classification, perceived own observability, or explicit bounds for missing data—will cement confidence in this parameter.


**Disentangling observability from message content**

The mechanism argument relies on the premise that the Calendar arm captures generic gift or reminder effects, while the Bracelet and Ink arms strictly isolate public observability (Sections 3.2 and 4.2.2).

Readers may note that the signaling treatments bundle observability with distinct messaging and associations. The bracelet explicitly carries a community-health and deworming message, whereas the calendar contains zero deworming message. Similarly, the ink mark is deeply familiar from voting contexts and entails private disutility. Consequently, differences between the Bracelet and Calendar arms could reflect variations in message content, campaign credibility, endorsement effects, or moral suasion across distances. Narrowing the causal claim to a public-signaling "bundle," or deploying belief, trust, and recall tests to isolate observability from message content, will clarify exactly what the experiment identifies.


**Local prevalence and the private valuation restriction**

In Section 6, the model ties the private utility of the calendar and bracelet together via a single, distance-invariant valuation parameter derived from Control communities, imposing $\beta_{Calendar}=\beta_{Bracelet}+\rho\psi$.

However, the descriptive evidence in Section 4.2.2 and Table B24 shows that gift choice among non-dewormers depends heavily on local item prevalence. When calendars are common, respondents choose bracelets at higher rates; when bracelets are common, the reverse is true. Because the structural valuation submodel does not incorporate this prevalence dependence, there is a distinct risk that the model attributes take-up differences to social image when they could mechanically stem from local saturation or scarcity overriding the fixed private-value restriction.


**Policy-invariant observability in the siting exercise**

The site-allocation counterfactual in Section 7 forms the basis for the paper's policy outcomes—changing the required sites from 108 to 81, with 19 sites attributed directly to the multiplier.

Applying the estimated demand function $T(z,d;\theta)$ directly to the planner's assignment introduces two strong assumptions. First, it requires extrapolating distance behavior up to 3.5 km, resting outside the experimental support limit of 2.5 km. Second, it assumes that the Control or Bracelet observability schedules depend solely on distance and signal regime, remaining fully invariant to the actual consolidation of sites, pooling of communities, and altered traffic patterns resulting from reducing the number of PoTs. Presenting sensitivity analyses that incorporate more pessimistic observability schedules under facility consolidation, and enforcing distance-cost caps strictly within experimental support, will ensure the site-count reductions carry the appropriate weight.

**Status**: [Pending]

---

## Detailed Comments (26)

### 1. Severe baseline imbalance in past deworming aligns with main effects

**Status**: [Pending]

**Quote**:
> The baseline sample is broadly balanced across item/signal arms. Statistically significant differences are sparse relative to the number of reported comparisons and do not reveal a systematic pattern across treatment arms. In particular, the social image measures in Panel C are balanced across arms.

**Feedback**:
Table C1 appears to undercut the statement that baseline differences do not reveal a systematic pattern. The variable "Dewormed in the past year" shows a large and statistically significant imbalance: Control is higher than Bracelet/Ink in Close communities but lower in Far communities, so the implied Far-Close differences for Bracelet and Ink are much less negative than in Control. Because prior deworming is plausibly predictive of current take-up, this imbalance should be addressed when interpreting the distance-by-public-signal take-up interactions, even if it does not by itself invalidate the randomized design.

---

### 2. Mathematically impossible sample means in balance tables

**Status**: [Pending]

**Quote**:
> | Age | 39.8 | 38.126 (0.758) | -0.714 [0.427] | -0.919 [0.286] | 0.232 [0.813] | 1.151 [0.208] | 0.5221 | 37.995 (0.838) | -0.804 [0.403] | -0.861 [0.324] | -0.798 [0.419] | 0.063 [0.945] | 0.7592 | 0.7776 |

**Feedback**:
In Table B1, the 'Sample mean' values for several variables fall completely outside the range implied by the 'Control mean' and treatment difference columns. For instance, the Sample mean for Age is 39.8, but combining the Control means with the corresponding treatment differences suggests all subgroup means correctly aggregate to a range between 37.1 and 38.4. The same pattern appears for the 'Female' and 'Phone owner' variables, and casades across Tables B2 and B12. This discrepancy apparent mechanically arises because the 'Control mean' column reports the intercept from the balance regression, representing the mean for the omitted county stratum rather than the raw overall Control mean. Presenting this omitted-stratum intercept as the generic 'Control mean' adjacent to the overall raw 'Sample mean' creates an apparent contradiction that is confusing for readers trying to deduce baseline characteristics from the table. To resolve this, consider either reporting raw Control means, employing fixed-effect adjusted means centered at the sample average of the fixed effects, or explicitly noting in the table notes that the 'Control mean' corresponds to the regression intercept for the omitted stratum.

---

### 3. Post-treatment selection in gift-choice test

**Status**: [Pending]

**Quote**:
> One might also worry that private value depends on local item prevalence. We therefore examine gift choices among non-dewormers in Bracelet and Calendar communities who did not already receive an item through deworming. Table B24 shows item prevalence effects: when the calendar is locally common, non-dewormers are more likely to choose the bracelet than in Control communities ( $+17.3 \mathrm{pp}, p<0.01$ ). When the bracelet is locally common, non-dewormers are less likely to choose the bracelet than in Control communities $(-9.7 \mathrm{pp}, p<0.05)$.

If marginal adults in Far Bracelet communities had a higher private valuation for the bracelet, bracelet choice should increase in Far Bracelet communities relative to comparable adults in Far Control communities. It does not: relative to Control-Far, bracelet choice in Bracelet-Far is lower by 8.9 percentage points ( $p<0.10$; Table B24). The Far-Close difference in the Bracelet coefficient is also small and statistically insignificant (+1.6 pp, s.e. 6.3 pp ).

**Feedback**:
In Section 4.2.2, the item-prevalence and Far Bracelet comparisons condition on being a non-dewormer, which is itself affected by the item/signal assignment. Since higher-bracelet-value adults in Far Bracelet communities may be precisely the people induced to deworm and then excluded from this sample, the -8.9 pp estimate does not by itself rule out high bracelet valuation among the relevant marginal compliers.

---

### 4. Unit error in relative valuation scale and reported estimates

**Status**: [Pending]

**Quote**:
> Let

$$
g_{i}^{\mathrm{val}} \equiv U_{c, i}-U_{b, i}
$$

denote respondent $i$ 's latent private value of the calendar relative to the bracelet, measured in the monetary units used in the elicitation. We assume

$$
g_{i}^{\mathrm{val}} \sim N(\psi, 1),
$$

where $\psi$ is the mean calendar minus bracelet value and the variance is normalized for scale.

**Feedback**:
There appears to be a unit inconsistency in the relative valuation model/reporting. Appendix J says $g_i^{\mathrm{val}}$ is measured in the monetary units used in the elicitation and fixes $g_i^{\mathrm{val}} \sim N(\psi,1)$, while Table B23 reports $\psi=0.472$ as KSh and also reports nondegenerate choice probabilities for 50 KSh offers. Furthermore, the text in Section 4.2.2 explicitly interprets the estimated mean gap as US$0.47, which implies a value about 100 times larger than 0.472 KSh. If $m_i$ entered the likelihood as raw KSh, a 50 KSh offer would be a 50-standard-deviation shift under the stated normalization; the reported probabilities instead suggest that the offer amounts were scaled, such as in dollars or hundreds of KSh. It is not fully clear whether $\psi$ is a monetary valuation or a normalized probit-index parameter, which matters for interpreting the valuation estimates and the mapping through $\rho$. The monetary scale used for $m_i$ and $\psi$ needs to be clarified.

---

### 5. Appendix E potential-distance check is underspecified

**Status**: [Pending]

**Quote**:
> A potential concern is that the acceptance-rejection algorithm may have introduced systematic differences in the set of feasible distances for communities assigned to Close versus Far, complicating the interpretation of realized continuous distance. Following Borusyak and Hull (2023), we assess this by simulating the realized selection and assignment mechanism. Holding the set of surveyed communities fixed, we re-run the full PoT selection and assignment algorithm 500 times and compute the centroid-to-PoT distances that would be induced under these counterfactual PoT-selection and assignment re-draws. We then compare the resulting potential distance distributions for communities assigned to Close and Far in the realized experiment. Figure E2 shows that these distributions are statistically indistinguishable, providing no evidence that the constrained algorithm generated systematic differences in potential distance profiles across realized assignment groups.

**Feedback**:
Appendix E’s potential-distance diagnostic is difficult to reconstruct from the description. In particular, it does not specify how each fixed surveyed community is assigned a counterfactual PoT distance when a simulated re-run selects a different set of PoT centers or targeted-community pairings, nor does it state the unit or weighting used for the Kolmogorov-Smirnov comparison.

---

### 6. Missing identification source for noise parameter

**Status**: [Pending]

**Quote**:
> Because latent utility is normalized by $\operatorname{Var}\left(v_{i}\right)=1, \sigma_{u}$ is interpreted as idiosyncratic decision noise relative to intrinsic motivation. The joint likelihood pins down the implied social image return, $\lambda p_{\text {Observed }}(z, d) \Delta\left(w^{*}\right)$, and how it changes with distance. These are the objects needed for the unitless social multiplier.

Remaining parameters. Conditional on the relative valuation parameter $\psi$, the noise parameter $\sigma_{u}$, and the observability parameters $\beta_{z}^{O}$ and $\gamma_{z}^{O}$, the remaining take-up parameters are identified through the take-up likelihood.

**Feedback**:
There seems to be an issue with the identification account for $\sigma_u$: the subsection explains how $\sigma_u$ affects the type-inference term and the multiplier, but it does not identify the empirical variation or parametric restriction that separately disciplines it. Because the later priors discussion notes that $\sigma_u$ requires regularization, the current identification paragraph leaves unclear whether this parameter is primarily identified by the nonlinear normal specification, by curvature in take-up across distance/observability cells, or materially by the prior.

---

### 7. Conflation of policymaker predictions with actual \mu=0 take-up

**Status**: [Pending]

**Quote**:
> Finally, consider a counterfactual in which the policymaker ignores social image and therefore believes $\mu=0$, so predicted demand is governed only by the estimated nonsocial component of demand, including travel costs. A practical analogue is a pilot or demand exercise based on private cost variation-for example, individual travel vouchers known only to recipients.

Such a design can estimate how private costs affect take-up, but it holds fixed the public information environment: it does not reveal how a community-level change in access costs alters observability, community beliefs, or the reputational payoff to participation. Under this counterfactual, the coverage target cannot be reached: even assigning every community to its nearest feasible PoT, using 144 sites and an average distance of 0.64 km , predicted take-up is only 9.1 percent (Table B30, Panel B).

**Feedback**:
The no-social-image counterfactual seems to be mixed with the forecast a policymaker would obtain from a private cost-variation pilot. If private vouchers vary travel costs while leaving the public information environment in place, the observed intercept would generally absorb the baseline social-image return; such data identify the private cost slope but not the nonsocial level of demand absent social image. Thus the 9.1% row is best interpreted as a structural $\mu=0$ demand counterfactual, not as the prediction from a policymaker who simply ignores social image when extrapolating from a private-voucher pilot.

---

### 8. Section 3.1 continuous-distance exogeneity needs support

**Status**: [Pending]

**Quote**:
> Distance as a continuous measure. Our primary reduced form distance treatment is the Close/Far assignment. We use realized community-centroid-to-PoT distance only in reduced form robustness checks, where it provides a continuous measure of the travel distance variation induced by the design. In the structural model, centroid-to-PoT distance is the main distance input, since it is the randomized spatial margin and a salient component of access costs. Specifically, we calculate the distance from each community's centroid to its assigned finalized PoT. All realized distance measures use the finalized PoT location at which treatment was delivered.

Because PoTs are selected using a geographically constrained procedure, the induced assignment mechanism could in principle depend on local geography, such as school density. Following Borusyak and Hull (2023), we re-run the full selection and assignment algorithm 500 times to assess whether communities assigned to Close and Far in the realized experiment had similar potential distance distributions under counterfactual re-draws. These simulations provide no evidence of systematic differences in potential centroid-to-PoT distances across the realized Close and Far groups (Figure E2). We further verify that realized centroid-to-PoT distance is uncorrelated with baseline covariates (Figure E3). Taken together, these checks support interpreting realized community-centroid-to-PoT distance as a continuous measure of the as-good-as-random travel distance variation in-

[^5]duced by the design.

**Feedback**:
The discussion of realized centroid-to-PoT distance could more sharply distinguish the randomized Close/Far assignment from the exact continuous distance realized after community selection, finalized PoT locations, and re-randomizations. The simulations and covariate checks support the interpretation, but it remains somewhat unclear what conditioning set justifies treating exact realized distance as as-good-as-random in the structural model.

---

### 9. Contradiction on prosocial signaling

**Status**: [Pending]

**Quote**:
> Notes: Green silicone bracelet distributed by CHVs at the point of treatment. The Swahili text reads "Tibu minyoo: boresha afya ya jamii yako", translated as "Treat worms: improve the health of your community".

**Feedback**:
There is a tension between Section 2's interpretation that the primary image motive is health conscientiousness/personal responsibility rather than prosociality and the intervention materials, which explicitly frame deworming as improving community health. The main social-image mechanism may still hold, but the exact image being signaled may include civic or prosocial content as well as personal responsibility.

---

### 10. Impossible crossing of curves in Figure 3(b)

**Status**: [Pending]

**Quote**:
> Panel (b) allows observability to fall with distance, $\mu^{\prime}(d)<0$. If the decline in $\mu(d)$ dominates the rise in $\Delta$, then the social image return can fall with distance: at large $d$, only highly motivated individuals deworm, but their action is less visible, so the reputational return $\mu(d) \Delta\left(w^{*}(d)\right)$ is small.

![](/documents/b1b14a42-0216-4727-8e8b-e1392684c111/images/image_002.jpg)
Figure 3: Equilibrium Social Image Returns $\mu(d) \Delta\left(w^{*}\right)$

**Feedback**:
Figure 3(b) appears to show the social-image-return curves crossing so that, for some positive $w^*$, a smaller idiosyncratic-noise dispersion such as $\sigma_u=0.5$ yields a lower $\mu(d)\Delta(w^*)$ than a larger dispersion such as $\sigma_u=1.5$. Under the stated normal model and a common declining observability schedule, $\Delta(w^*)$ should fall as $\sigma_u$ rises at any fixed cutoff, and the implicit equilibrium mapping should preserve this ordering when $\mu'(d)\le 0$. The lower-$\sigma_u$ curve should therefore remain weakly above the higher-$\sigma_u$ curve; the apparent crossing suggests a possible inconsistency in how Panel (b) maps cutoff, distance, and observability.

---

### 11. Outdated $\hat{R}$ convergence threshold for HMC

**Status**: [Pending]

**Quote**:
> The model is estimated in Stan using Hamiltonian Monte Carlo with four chains, 400 warm up draws, and 400 saved draws per chain. Sampler diagnostics show no divergent transitions and split $\hat{R}<1.1$ for all parameters.

**Feedback**:
The HMC diagnostics are reported using a relatively loose convergence threshold. For a Bayesian structural model that feeds into decomposition and policy counterfactuals, $\hat{R}<1.1$ alone does not fully rule out mixing problems, especially with 400 saved draws per chain. If the actual values are all comfortably below stricter thresholds and effective sample sizes are adequate, this is only a reporting issue; otherwise, the posterior summaries may require longer runs or additional convergence evidence.

---

### 12. Appendix E.1 permutation test may not match the design

**Status**: [Pending]

**Quote**:
> As an additional check that our continuous distance measure is orthogonal to pretreatment characteristics, we regress each baseline covariate on the distance from the community centroid to its assigned point of treatment (PoT), controlling for county fixed effects and clustering standard errors at the community level. We then permute community-level distances across communities within county strata 500 times and recompute the associated test statistic. Figure E3 plots the resulting randomization distribution for each covariate under the null of no association. The vertical line marks the realized value, and permutation $p$-values based on the 500 draws are reported in the upper right-hand corner of each panel. For all covariates shown, the $p$-value exceeds 0.10 , indicating that realized centroid-to-PoT distance is not systematically related to observable characteristics.

**Feedback**:
Appendix E.1’s permutation exercise is useful as a descriptive balance check, but it is not clear that freely permuting realized community-level distances within county reproduces the constrained PoT-selection and community-assignment mechanism. If the reported values are intended as design-based randomization-inference $p$-values, the exchangeability assumption behind this permutation scheme needs to be justified.

---

### 13. Section 3.2.1 SMS counts are hard to reconcile

**Status**: [Pending]

**Quote**:
> Specifically, we randomly sampled 17 phone-owning adults for Reminder Only and 17 for Social Information in Control clusters, and 29 phone-owning adults for Social Information in clusters assigned to Bracelet, Ink, or Calendar. Where possible, we also identified within-household backup phone owners to use when the initially sampled adult could not be reached. Only one adult per household could ultimately be enrolled in an SMS treatment. Enumerators recruited sampled individuals in person prior to the campaign. Only reached and consenting adults could receive SMS messages. ${ }^{12}$

**Feedback**:
There seems to be an issue with the SMS sample counts: applying the stated per-cluster sampling rule to 34 Control clusters and 110 non-Control clusters implies 4,346 primary SMS sampling slots, while later descriptions report 3,022 enrolled SMS recipients and the appendix appears to use 3,743 initially sampled or assigned SMS adults. The text does not clearly distinguish target slots, actual primary draws, backups, reached individuals, consenting/enrolled recipients, and monitored SMS adults.

---

### 14. Distance-related omission needs attention in Section 3.5

**Status**: [Pending]

**Quote**:
> Table B3 examines whether item/signal assignment predicts omission from the monitored take-up sample. The dependent variable is an indicator equal to one if an individual assigned to monitoring was omitted from the PoT monitoring list. Across the Bracelet, Calendar, and Ink arms, the estimated coefficients are small and statistically indistinguishable from zero, both in the pooled specification and when the sample is split by distance. These results provide no evidence that item/signal assignment induced differential omission from the monitored take-up sample.

**Feedback**:
In Section 3.5, the text notes that item/signal assignment does not predict omission from the administrative monitoring list. However, it does not discuss omission by distance, which is a core treatment dimension. As the reviewer rightly notes, Appendix Table B3 reveals a statistically significant difference in the Control row of the Far-Close column (-0.014, $p < 0.05$), meaning omission is about 1.4 percentage points lower in Far Control communities than in Close Control communities. While this 1.4 pp magnitude is likely too small to drive the large main take-up effects of distance (e.g., the ~16.2 pp gap), leaving a statistically significant imbalance in attrition unaddressed in the text can invite unnecessary sample-selection concerns. Adding a brief explanation for this distance-related imbalance (and perhaps noting why its small magnitude or mechanical nature renders it harmless) would fully close this loop for the reader, similar to the careful discussion already provided for the observability module missingness.

---

### 15. Confusing description of re-randomization

**Status**: [Pending]

**Quote**:
> Eligibility is then verified using the distance between the selected community's centroid and the finalized PoT. If the finalized PoT location or initially selected community made the assigned distance band infeasible, the PoT-community pairing was re-randomized before implementation using the same protocol and eligibility criteria.

**Feedback**:
The description of the six pre-implementation re-randomizations is somewhat ambiguous. The realized counts move from 74 Close/70 Far after attrition to 80 Close/64 Far, and Appendix E later says all six re-randomized cases were realized as Close assignments; this suggests that the Close/Far assignment, not merely the community drawn for a fixed assignment, changed in those cases. It would be useful to make clear what part of the PoT-community-distance assignment was re-randomized and how feasibility restrictions entered that re-randomization.

---

### 16. Inconsistent observability estimates in Introduction

**Status**: [Pending]

**Quote**:
> We first document how distance and the public signal treatments affect observability. In the no-item Control arm, deworming is highly observable in Close communities: respondents report knowing 75 percent of peers' deworming status. Observability declines with distance: reported knowledge is 12 percentage points lower in Far than in Close Control communities ( $p=0.056$ ), consistent with learning through direct encounters at treatment sites or seeing others en route to treatment. The public signals largely offset this distance related decline. In Far communities, relative to Control, bracelets raise knowledge of peers' status by 22 percentage points ( $p<0.001$ ) and ink raises it by 16 percentage points ( $p=0.007$ ), while in Close communities, the same signals raise observability by only $3-5$ percentage points ( $p=0.357-0.644$ ).

**Feedback**:
The Introduction’s observability summary appears not to match the main estimates later in Section 4.1.1/Table 1. For example, Table 1 reports a Control Far-Close observability difference of 14.6 percentage points, Bracelet and Ink effects in Far communities of 24.8 and 16.9 percentage points, and Close-community effects of 6.4 and 5.2 percentage points. The qualitative message is unchanged, but the introductory point estimates and reported $p$-values should be checked against the final specification.

---

### 17. Contradiction between text and Figure M1

**Status**: [Pending]

**Quote**:
> ![](/documents/b1b14a42-0216-4727-8e8b-e1392684c111/images/image_021.jpg)
Figure M1: Simulation under Gaussian and Bimodal Type Distributions

Notes: Panel (a) plots the assumed densities for $w$: Gaussian in red and bimodal in blue.

**Feedback**:
There appears to be a disconnect between the simulation parameters described in the text of Section M and those plotted in Figure M1. The text specifies a standard normal $w \sim N(0,1)$ and an "equally weighted bimodal mixture of $N(-\eta, 1)$ and $N(\eta, 1)$." However, Panel (a) of Figure M1 displays a highly asymmetric bimodal distribution (modes at approximately 0.4 and -2.0, with differing peak heights) and a Gaussian benchmark with a peak density of 2.0, which implies a standard deviation of approximately 0.2 rather than 1. Additionally, the axes in panels (a) and (c) are labeled $v^*$ and $V^*$, rather than $w$ or $w^*$ as denoted in the text and caption. Updating the text to accurately reflect the simulation parameters plotted in the figure would resolve this confusion.

---

### 18. Apparent contradiction regarding closest PoT

**Status**: [Pending]

**Quote**:
> Study communities are selected from the remaining, non-overlapping portion of the 2.5 km catchment. As a result, any selected study community is at least 0.5 km closer to its assigned PoT than to the nearest alternative PoT.

**Feedback**:
The 0.5 km separation statement appears to refer to the candidate PoT centers used in the spatial selection algorithm, but the later statement that in two clusters the assigned finalized PoT was not the closest deworming site creates a tension as written. The distinction between candidate PoT centers and finalized PoT locations should be made explicit at the point where the 0.5 km guarantee is stated.

---

### 19. Sign error in reported Far-Close difference

**Status**: [Pending]

**Quote**:
> Control respondents report knowing the status of 68.1 percent of recognized peers (Table 1). This share is 74.5 percent in Close communities and 59.9 percent in Far communities, a Far-Close difference of 14.6 percentage points.

**Feedback**:
The reported Far-Close difference appears to have the wrong sign: 59.9 percent in Far communities minus 74.5 percent in Close communities is -14.6 percentage points, consistent with the negative Far-Close estimate in Table 1.

---

### 20. A6 outreach-script note appears to describe A7 flyers

**Status**: [Pending]

**Quote**:
> 7. You will receive $\_\_\_\_$ for deworming yourself as a symbol of your passion towards improving the health of the members in your family and the community.

## Figure A6: Community Outreach Script

Notes: This figure displays the four one-page flyers that CHVs distributed during the pre-campaign information visits. Each flyer, printed in Swahili, announces the free adult deworming campaign, leaves blank fields for the locally customized site ("Wapi") and date ("Lini"), and, except in the Control arm, depicts the item or signal offered upon deworming: bracelet, calendar, or ink mark. Apart from the item/signal panel, wording and layout are identical across flyers.

![](/documents/b1b14a42-0216-4727-8e8b-e1392684c111/images/image_011.jpg)
Figure A7: Flyers

Notes: This figure displays the four one-page flyers that CHVs handed out during the pre-deworming information visits. Each flyer-printed in Swahili-announces the free adult deworming campaign, leaves blank fields for the locally customized site ("Wapi") and date ("Lini"), and, except in the Control arm, depicts the incentive offered upon deworming (bracelet, calendar, or ink mark). Apart from the incentive panel, wording and layout are identical across flyers to preserve comparability.

**Feedback**:
The note under Figure A6 appears to be misassigned: it says the figure displays four flyers with “Wapi”/“Lini” fields and item panels, which describes Figure A7, while Figure A6 is titled and preceded by the community outreach script. This creates a factual mismatch in the appendix documentation.

---

### 21. Missing conclusion about calendar and bracelet

**Status**: [Pending]

**Quote**:
> The parameter $\psi$ helps pin down the relative private value of the calendar and bracelet. The
estimates imply a positive weight on social image, a positive distance cost, and a negative private value for ink.

**Feedback**:
In Section 6.1, the paragraph introduces $\psi$ to capture the relative private value of the calendar and bracelet. However, the subsequent summary sentence omits the empirical finding for $\psi$ and instead states a finding regarding ink's private value, which was not among the previously listed parameters. Briefly stating the estimate for $\psi$ and aligning the summary with the listed parameters would resolve this local disconnect.

---

### 22. Arm variation is not represented in the Section 5 model

**Status**: [Pending]

**Quote**:
> Utility from action $y_{i}$ is given by

$$
U_{i}\left(y_{i}\right)=\left(b-d+v_{i}+u_{i}\right) y_{i}+\mu(d) E_{-i}\left[V \mid y_{i}\right]
$$

where $b-d$ is the net private payoff, combining the health and material item value $b$ with the distance induced travel cost $d \geq 0$. Intrinsic motivation $v_{i}$, interpreted as health conscientiousness or civic mindedness, is privately known and unobservable to others. The idiosyncratic shock $u_{i}$ captures taste variation or one time costs and is independent of $v_{i}$. Since both variables are latent, we define the composite type $w_{i} \equiv v_{i}+u_{i}$.

The term $\mu(d) E_{-i}\left[V \mid y_{i}\right]$ captures social image. Following Benabou and Tirole (2025), $\mu(d)$ bundles the physical observability of the act and the individual's valuation of social image. In our context, observability is empirically lower in Far than in Close communities absent public signals, so we allow $\mu(d)$ to vary with distance and with the item/signal arm.

**Feedback**:
There seems to be an issue with the notation in Section 5.1: the displayed utility uses $\mu(d)$ and a single private payoff $b$, even though the surrounding text says observability and item/signal arms can vary by treatment arm. The later structural specification resolves this with arm-specific private-payoff and observability terms, so the concern is local, but the introductory formal model suppresses treatment-arm dependence in a way that obscures the mapping from the reduced-form arms to the model primitives.

---

### 23. Appendix F ambiguous use of “controls”

**Status**: [Pending]

**Quote**:
> Data collection issues affecting analysis samples. Two programming or list generation issues affect the construction of analysis samples. First, 255 respondents assigned to the observability module did not receive it because of a programming error. Section 3.5 discusses this issue, Table B4 reports differential missingness checks, and Table B6 reports Lee bounds for the observability outcomes. Second, a SurveyCTO list generation error caused 989 wave 1 non-phone owner controls to be omitted from the PoT monitoring lists. These individuals are therefore absent from the administrative take-up data. Section 3.5 discusses this issue, Table B3 tests for differential omission across item/signal arms, and Appendix $D$ describes the construction of the administrative monitoring sample.

**Feedback**:
There seems to be an issue with the term “controls” in Appendix F: in this sentence it appears to refer to no-SMS comparison controls rather than the no-item Control arm, but the paper also uses “Control” for the no-item arm. Because the sentence concerns monitoring-list omissions, this ambiguity can affect how readers interpret the sample-construction problem.

---

### 24. Incorrect statistical terminology: normal mixture

**Status**: [Pending]

**Quote**:
> We normalize $v_{i} \sim N(0,1)$ and assume $u_{i} \sim N\left(0, \sigma_{u}^{2}\right)$, so $w_{i}=v_{i}+u_{i} \sim N\left(0,1+\sigma_{u}^{2}\right)$. For the normal mixture, the type inference term is

**Feedback**:
The phrase “normal mixture” is statistically imprecise here. Since $w_i=v_i+u_i$ is defined as the sum of independent normal variables and is immediately stated to be normal, this is an additive-normal/convolution specification rather than a mixture distribution.

---

### 25. Possible error in standard error reporting in Table B21

**Status**: [Pending]

**Quote**:
> | Social Info | $0.131^{* * *}$ |
|  | $(0.028)$ |
| Control mean | 0.330 |
|  | $(0.470)$ |
| Observations | 1,705 |


Notes: This table reports estimates from a linear probability model (OLS) with a binary outcome for deworming take-up. The regressors are indicators for receiving a Social Information SMS or a Reminder Only SMS. The omitted group is phone owners who received no SMS. Because the Reminder treatment was implemented only in Control communities, the analysis sample is restricted to phone owners in the Control arm, not the full SMS treatment sample ( $N=3,022$ ). All regressions include county fixed effects and the LASSO selected controls for sex and age, and expected distance as a design control. Bootstrapped standard errors, clustered at the community level, are in parentheses. ${ }^{*} p<0.10,{ }^{* *} p<0.05,{ }^{* * *} p<0.01$.

**Feedback**:
In Table B21, the parenthetical value below the Control mean appears to be the control-group standard deviation of the binary take-up outcome, not a bootstrapped clustered standard error. This conflicts with the note stating that values in parentheses are standard errors.

---

### 26. Discrepancy between text and Table C1 statistics

**Status**: [Pending]

**Quote**:
> When asked who is at risk of worm infection, 94 percent name children and 70 percent name adults.

Adult take-up and program context. Despite low monetary costs-a full course costs roughly US\$0.50-2.00 in local pharmacies-regular adult deworming remains limited: 38 percent report deworming in the previous year, although 69 percent have done so at least once (Table C1, Panel B).

**Feedback**:
The whole-percentage summaries in Section 2 appear slightly inconsistent with Table C1, Panel B. Table C1 reports 0.934 for knowing children get worms, 0.374 for dewormed in the past year, and 0.684 for ever dewormed, which would conventionally round to 93%, 37%, and 68%, rather than 94%, 38%, and 69%. The differences are small, but the prose and appendix should be consistent.

---
