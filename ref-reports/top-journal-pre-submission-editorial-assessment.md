# Overall assessment

This is a **serious, original, and unusually ambitious paper**. The project has a top-journal question, a large and creative field experiment, administrative outcome data, and an elegant conceptual mechanism. The draft is substantially stronger than the typical working paper.

My editorial summary would be:

> **Top-journal-caliber idea and data; not yet top-journal-caliber identification of the headline parameter.**

The reduced-form paper is compelling. The precise claims that the social multiplier equals roughly 0.8–0.9 or 1.2–1.6, and that it changes the required number of treatment sites from 108 to 81, are materially more assumption-dependent than the abstract and introduction currently convey.

At a top-five general-interest journal, I would regard this as a **borderline but plausible send-out**, not an obvious desk rejection. My modal forecast after refereeing, however, would currently be rejection rather than an R&R. At a leading field journal, it would already be a strong R&R candidate.

| Dimension | Editorial assessment |
|---|---|
| Economic question | Excellent |
| Experimental scale and ambition | Excellent |
| Reduced-form causal evidence | Strong, with one important weak point |
| Social-image interpretation | Plausible and well triangulated, but not decisive |
| Structural quantification | Ambitious; currently the weakest part |
| Policy application | Interesting illustration, but not yet decision-grade |
| Writing and organization | Clear and professional, but too long and not fully finished |
| Current top-five readiness | Promising but below the bar |

## What is particularly strong

### 1. The economic mechanism is genuinely interesting

The central idea is broader than deworming: lowering an access cost changes not only participation but also what participation signals about the participant. That is a clean, portable economic mechanism. It connects behavioral public economics, development, health, social interactions, and spatial policy.

The distinction between two forces is especially appealing:

- Greater distance can make participation a stronger signal of motivation.
- Greater distance can simultaneously make participation less observable.

That these forces can produce either amplification or mitigation gives the paper a clear conceptual contribution rather than merely another intervention-effects result.

### 2. The experiment is unusually well designed and implemented

The cross-randomization of distance and item/signal assignment is clever. The paper also has several features editors value:

- 144 randomized communities and almost 10,000 administratively monitored adults.
- An objectively measured primary outcome rather than self-reported take-up.
- Two public signals with different properties.
- A tangible comparison item rather than only a no-treatment control.
- Auxiliary SMS interventions.
- Measures of reported peer knowledge, perceived observability, gift choices, and relative valuations.
- Extensive balance, attrition, implementation, and robustness analyses.

The large distance effect in the control arm—about a 16 percentage-point decline—is clean and economically substantial. The effects on reported observability are also quite persuasive: public signals matter much more in Far communities, exactly where background observability is lower. The bracelet-versus-calendar observability comparisons are particularly useful because they move the analysis beyond a simple gift effect.

### 3. The appendices reflect serious empirical craftsmanship

The draft is transparent about post-PAP additions, the distance reclassification issue, the monitoring-list error, missing observability modules, and sample construction. The ITT specifications, Lee bounds, alternative observability measures, continuous-distance checks, and robustness models are all valuable.

An editor would notice that the authors have anticipated many obvious objections rather than hiding them. That materially raises the paper’s credibility.

---

# The issues most likely to determine the editorial decision

## 1. The central reduced-form take-up interaction is less decisive than the framing suggests

The cleanest causal findings are:

1. Distance sharply reduces take-up.
2. Public signals sharply increase reported observability in Far communities.
3. Bracelets increase take-up, especially in Far communities.

The evidence that public signals **attenuate the distance gradient in take-up**, however, is more modest.

In the primary binary Close/Far specification, the pooled Signal-versus-No-Signal Far–Close difference is 6.8 percentage points with \(p=0.083\). The bracelet-calendar difference in the interaction has \(p=0.101\). Under continuous distance, the pooled interaction reaches \(p=0.049\), but the bracelet-calendar interaction has \(p=0.336\). The original-assignment ITT results are stronger, but that makes the choice of the principal specification consequential.

This is not a fatal weakness. The signs are generally consistent, and the observability results are strong. But it affects how aggressively the paper can describe the multiplier as an experimental finding. A skeptical editor will ask:

> Is the paper’s main result a statistically decisive treatment interaction, or is it a suggestive interaction that is made precise by the structural model?

At present, the latter is closer to the truth.

I would make the original randomized assignment—or an IV specification using original assignment to instrument realized distance—the principal design-based analysis. The paper should also report randomization inference or wild-cluster-bootstrap inference for the main omnibus tests and be very explicit about which pooling and distance specifications were preregistered.

## 2. The experiment changes more than observability

The bracelet, ink, calendar, and control arms are not pure manipulations of the probability that a fixed action is observed.

The bracelet:

- Has private consumption and identity value.
- Carries an explicit community-health message.
- May work as a reminder or commitment device.
- May alter the meaning of participation, not simply its visibility.

The ink mark has its own private disutility and voting-related familiarity. The calendar is widely retained and displayed at home and is often understood as being connected to deworming, even though it is less publicly portable.

The calendar comparison, gift choices, WTA exercise, and SMS interventions all help. But they do not fully recreate the clean counterfactual of **the same signal being privately versus publicly observable**. Nor do the SMS treatments conclusively eliminate descriptive norms, social learning, or social proof: seeing a neighbor with a bracelet is not necessarily behaviorally equivalent to receiving an aggregate take-up statistic by text.

The strongest interpretation supported by the design is therefore:

> Publicly legible participation signals interact with access costs in a manner consistent with social-image incentives.

That is somewhat weaker than:

> The experiment cleanly identifies the social-image mechanism and its equilibrium multiplier.

Without a new public-versus-private version of the same item, this limitation may not be completely fixable. The right response is either to temper the mechanism claim or to provide formal bounds that permit signal-specific private utility and signal-specific reputational content.

## 3. The observability measure does not exactly match the theoretical object

The main empirical observability outcome is the share of recognized peers for whom a respondent gives a definite Yes or No answer. This is an informative measure, but it is not literally the probability that an actor expects their action to be observed.

There are three related concerns.

First, conditional accuracy falls in some signal cells as “Don’t know” answers are converted into definite classifications. In Table B16, signals increase unconditional correct classification in Far communities, which is reassuring, but the added classifications are not all genuine observation. Some may be inference or guessing.

Second, the theoretically closest measure is respondents’ belief that others know **their own** deworming status. Those estimates are directionally supportive but appreciably less precise than the primary peer-knowledge measure in Table B17.

Third, bracelets and ink are essentially one-sided, noisy signals of participation. The structural model treats the action as observed with probability \(p_{\text{Observed}}\) and then applies the same image structure to participation and abstention. In reality:

- A visible bracelet is strong evidence of participation.
- The absence of a bracelet is not necessarily evidence of abstention.
- Wearing and retention are incomplete.
- Ink fades.
- False positives and false negatives need not be symmetric.

A top-journal referee is likely to request an explicit noisy-signal model with separate observation technologies conditional on participation and nonparticipation. At minimum, the paper should show how the multiplier changes using perceived observability and correct classification rather than definite-answer rates.

## 4. The structural model is doing more work than the paper acknowledges

Section 6 candidly states that the pure travel-cost effect is not directly observed: every experimental arm has positive observability. The multiplier is consequently recovered through cross-arm restrictions—most importantly a common distance-cost parameter and a common social-image weight—combined with arm-specific observability schedules and private-payoff restrictions.

That is legitimate structural economics. But it creates four major editorial concerns.

### The common-image-weight restriction is strong

The model assumes signals change the probability of observation but not the reputational value of being observed, conditional on the action. Yet the bracelet explicitly communicates community-mindedness, while ink and direct observation at a treatment site may carry different meanings.

If the bracelet changes both \(p_{\text{Observed}}\) and the image weight \(\lambda\), the cross-arm distance gradients no longer cleanly identify the decomposition. The paper needs a serious sensitivity analysis allowing \(\lambda\) to differ across signals, or partial-identification results showing what can be learned under bounded heterogeneity.

### The valuation data barely discipline the multiplier

The prose assigns an important identification role to the bracelet-calendar valuation exercise. But Table B28 indicates that the WTP likelihood has almost no local influence on the social multiplier, and the estimated mapping parameter \(\rho\) is practically near zero. In effect, the auxiliary valuation data and the take-up utility scale are only weakly connected.

That makes the claim that the valuation exercise “bounds” the private-value explanation stronger than the estimated model itself appears to support.

### The uncertainty calculation appears too optimistic

Treatment varies at 144 communities. The structural likelihood, as presented, uses individual Bernoulli take-up observations and respondent-level binomial observability observations. I do not see a cluster-level disturbance, a cluster bootstrap of the full estimation, or a design-based correction for within-community dependence.

This matters because the structural model produces much tighter inference than the reduced form. For example, the reduced-form pooled interaction is marginal, while the structural counterpart has a narrow credible interval comfortably excluding zero. Some tightening is expected from model restrictions, but the magnitude of the gain raises a concern about pseudo-replication.

A top-journal version should re-estimate the model through community-level bootstrap resampling, introduce hierarchical community shocks, or use a minimum-distance procedure based on cluster-level experimental moments and their clustered covariance matrix. Overlap among the take-up, observability, and valuation samples should also be reflected in uncertainty.

### The headline depends materially on the assumed information structure

The full-information specification does not merely make the multiplier somewhat smaller. In Table L5 it eliminates, and in places reverses, the amplification result for the control environment. The paper discounts this model because it does not fit the Far–Close treatment patterns as well. That is relevant, but not sufficient to establish that observers use the community-centroid distance assumed by the preferred model.

The ideal evidence would directly show what peers know about one another’s access costs. In its absence, the paper should treat the information structure as a source of identification uncertainty rather than simply selecting the model that fits best.

The broader editorial point is:

> A reader can believe all the randomized reduced-form evidence and still remain unconvinced by the precise 0.8–1.6 multiplier estimates.

That gap is the manuscript’s central vulnerability.

## 5. The policy exercise currently outruns the evidence

The site-allocation exercise is clever and makes the multiplier tangible. But it is better viewed as an illustrative model application than as a reliable policy prescription.

The exercise:

- Extrapolates beyond the experimental support to a 3.5 km cap.
- Minimizes the number of sites rather than total social cost or welfare.
- Omits household travel costs from the objective.
- Does not price the signal intervention, training, supervision, or crowding.
- Uses a common coverage target and equal community weights.
- Assumes the estimated observability schedule remains valid when several communities are pooled at fewer sites.
- Assumes the meaning and visibility of bracelets remain stable at a different program scale and site configuration.

The last point is especially important. Moving from 144 sites to 81 changes the pattern of social encounters and which communities mix at a treatment site. The observability function estimated under the experimental allocation need not transport unchanged to the optimized allocation.

The 108-versus-81 result is memorable, but because it rests on both the structural multiplier and several allocation assumptions, it should not carry as much weight in the abstract as it currently does. Table B35 is useful; it does not yet constitute a policy-ready counterfactual.

---

# Presentation and manuscript readiness

The paper is generally clear, but it is too long and repetitive for a top-journal submission. The introduction repeatedly explains the same observability-versus-inference decomposition, and the main text could probably be shortened by 20–25 percent without losing substance.

There are also visible signs that this is not the final submission draft:

- Unresolved citation placeholders appear in the site-allocation discussion.
- Table F1 contains an unresolved “??”.
- Some numerical descriptions in the introduction do not exactly match the principal tables.
- The 80-versus-81 site presentation needs a transparent explanation.
- The notation for the private payoff coefficients in the structural model and Table B29 is not fully intuitive.
- “Observability,” “reported knowledge,” “correct classification,” and “perceived observability” should be kept more sharply distinct.

These are fixable, but at a top journal they matter because they make an already complicated identification argument harder to audit.

# What I would require before considering an R&R

1. **Put the randomization back at the center.** Make the original Close/Far assignment and a design-based IV analysis primary. Report cell means, cluster-level plots, randomization inference, wild-cluster-bootstrap results, and appropriate multiple-testing adjustments.

2. **Narrow or strengthen the mechanism claim.** Use perceived own-observability and correct-classification measures prominently; analyze campaign-day dynamics; explicitly address descriptive norms and social learning; and distinguish social image from other public-signal mechanisms. Without new data, describe the evidence as strongly consistent with social image rather than uniquely identifying it.

3. **Rebuild structural uncertainty and sensitivity.** Use cluster-resampled or minimum-distance estimation, allow signal-specific image weights, model noisy and asymmetric observation, provide broad prior sensitivity, and report profile or identified-set information for the multiplier.

4. **Treat the direct motivation-inference channel more cautiously.** It is currently inferred from the latent-type structure rather than directly measured. The draft should show how much of the result survives across estimated alternative type distributions and information structures, not only illustrative simulations.

5. **Either strengthen or demote the site-allocation exercise.** Include signal costs, household travel costs, population weighting, and alternative observability-under-pooling assumptions—or reposition the exercise as an appendix illustration rather than a headline policy result.

# Bottom line

This is a **high-quality paper with a real chance of becoming a top-journal paper**. Its strongest publishable core is the experimental result that access costs and public visibility interact: distance changes both participation and the informational environment surrounding participation.

The current manuscript’s weakness is not lack of ingenuity or data. It is that the paper turns a persuasive but moderately precise experimental pattern into a very precise structural multiplier and a striking spatial-policy result. The structural assumptions, observability measurement, and uncertainty treatment are not yet strong enough to bear that weight.

My editorial disposition would therefore be:

> **Borderline top-five send-out; likely rejection in its present form, but with unusually high upside. Strong R&R potential if the paper makes the structural claims more transparent, more design-based, and less precise than the data warrant.**
