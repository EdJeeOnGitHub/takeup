# Referee Concerns: Sample Definition and PAP Deviations


## 1. Sample Definition / Creation Confusion

### QJE R1

> "The randomization protocol and sampling frame are complex. It took me a few rounds of reading to understand the protocol and I am still not sure I fully grasped it. I suggest a clear schematic illustration of the randomization hierarchy and the relationship between different samples (take-up, baseline, observability, endline, etc.). It is currently difficult to track how these samples overlap or why, for example, the endline take-up sample is larger than the observability sample."

> "I find the description of the accept–reject randomization algorithm also confusing. I read the PAP and I got confused even more as I couldn't reconcile every aspect of the two expositions."

**Addressed:** We have added:

- Sample construction + flow diagram added. Explicit links between endline modules/samples discussed in appendix + main text.
- The accept–reject algorithm is now formally described in the appendix subsection with explicit notation.

**Not addressed:** 

### QJE R3

Under "Minor Comments":

> "I didn't understand the paragraph on cluster selection and assignment. What does 'drawing 158 clusters' mean exactly? What is the definition of a cluster here?"

> "What does this mean: 'For each cluster, we designated an "anchor school" that satisfied the assigned distance category'? I guess a cluster is either Close or Far, so if a cluster is Close, you are selecting an anchor school that is 0-1.25 km away from... what?"

> "Appendix C and Figure A3 helped a bit. But still, the main text could be clearer, I think. For example, the selection of 'clusters' could be more clearly described as the selection of 158 primary schools to be points of treatment, subject to these selected schools having non-overlapping X km radii."

> "Though then I don't understand what it means for clusters to be randomly assigned to close vs. far 'from the cluster's centroid to the PoT' (from the appendix). I had thought that the centroid *was* the PoT (the primary school)."

> "One detail I missed through not understanding: is a given community always assigned the nearest school as the PoT, with the close vs. far coming from choosing a PoT that is near the centroid of a cluster vs. near the edge? Or does the variation come from whether a community is randomized to a close vs. far PoT? I guess it must be the former because the randomization is at the PoT-level rather than the community-level?"

Also notes that it appears some data between 2.5 km and 3.5 km is dropped from some analyses (see Figures 3 and A3), raising a question about sample consistency.

**Addressed:** 
- Cluster definition (2.5 km catchment + one targeted community) is now stated explicitly.

- Use PoT as shorthand for cluster, apart from in appendix where full detail added.

- Remove anchor discussion and only mention in appendix - people don't need every operational detail here I think.

- "158 clusters" is explained as 150 targets plus 8 fallbacks, with implementation limited to 144. 
- The centroid/PoT confusion is resolved: the distance measure is the community centroid to the PoT, not PoT to PoT. 
- The anchor school is explained as a field reference point for enumerators locating nearby communities, not a structural design element. 
- Fig A2 distance cut-off: explained by kernel density smoothing, no one is dropped from the analysis

**Not addressed:** 


### AER Ref1

Raises essentially the same cluster/design confusion as QJE R3:

> "I didn't understand the paragraph on cluster selection and assignment. What does 'drawing 158 clusters' mean exactly? What is the definition of a cluster here?"

> "What does this mean: 'For each cluster, we designated an "anchor school" that satisfied the assigned distance category'?"

> "Though then I don't understand what it means for clusters to be randomly assigned to close vs. far 'from the cluster's centroid to the PoT' (from the appendix). I had thought that the centroid was the PoT (the primary school)."

> "One detail I missed through not understanding: is a given community always assigned the nearest school as the PoT... Or does the variation come from whether a community is randomized to a close vs. far PoT?"

Sample size discrepancy:

> "Baseline survey data: 144 communities * 15 households * 1 adult = 2,160 ≠ 4,823 adults? Where does the discrepancy in numbers come from?"

Also asks what the distance treatment *is* — whether PoTs were created for the study or pre-existing, and whether the randomization assigned close vs. far to communities or to facilities.

**Addressed:** 

- PoTs are pre-existing primary schools/religious centers (stated in main text and appendix definitions paragraph). 
- The centroid/PoT distinction and the 3-step explanation should make randomization mechanism clearer. 
- Baseline N sample discrepancy - sample flow diagram/discussion should cover.
- Anchor school moved to appendix.

**Not addressed:** 

### AER Ref2

> "There are several discrepancies between the project as described in the AEA Registry and in the paper. For instance, the number of treatment units is 150 in one case and 144 in another. Perhaps even more importantly, while the paper states that 'Communities are randomly assigned to treatment locations either close or farther away,' the registry does not mention distance as a randomized treatment arm. Instead, similar to the Sierra Leone study, this is described more as a sampling decision: 'One community per central location catchment area was randomly selected, with half of the communities being less than 1.25km from the central location and half of the communities being more than 1.25km away located from the central location.'"

> "The description of the main outcome measurement is somewhat unclear. Specifically, what does 'monitored at the point of treatment' mean in 'A sample of 9,805 adults whose deworming status was monitored at the point of treatment by enumerators'? It seems to suggest potential selection into measurement. Why would people show up at these points if not for the treatment?"

> "You do not mention the size of the treatment cells interacted with the distance heterogeneity, nor the sample size within each cell. Providing this information would add clarity."

**Addressed**:

- Explicit main text/appendix text about number of PoTs 158 -> 144.
- Explicitly discuss distance randomization.
- Treatment cell size now explicit.

**Not addressed:**

- Distance randomization not explicit in PAP?



### AER Ref3

> (same cluster/PoT/centroid confusion as QJE R3 and AER Ref1 — see those entries)

> "It seems that rather than close and far, the distance measure is continuous (see Figures 3 and A3 — side point: it appears that there is some data between 2.5km and 3.5km that is dropped for some of the analysis). The 1.25km cutoff seems arbitrary... Could the authors show what happens when this cutoff is changed to be bigger or smaller?"

**Addressed:** 
- The 1.25 km cutoff is justified as chosen ex ante — the midpoint of the 2.5 km catchment radius and approximately half the upper range of willingness-to-travel reported in piloting.

- The continuous distance measure and its use alongside the binary indicator are now explicitly described. 

- Data past 2.5km discussed.


**Not addressed:**  

- The request for cutoff sensitivity analysis is not addressed (but I don't think we'd want to re-define treatment status)

---

## 2. PAP / Deviations from Pre-Analysis Plan

### QJE R1

> "Much of my concern arises from the fact that the experiment was not designed ex ante for structural estimation of social multipliers, and the main treatments and conceptual framing diverge significantly from what was pre-registered in the pre-analysis plan (PAP). The paper's modeling exercise is interesting in its own right, but it feels somewhat ex post—an attempt to fit the data to the Benabou–Tirole framework rather than a context naturally suited to it."

> "I read the PAP and I got confused even more as I couldn't reconcile every aspect of the two expositions [of the randomization protocol]."

**Addressed:** 

- The appendix explicitly states "The notation used matches that used in the Cluster Selection Algorithm section 4.4 of the pre-analysis plan," which should help reconcile the two expositions.


- Structural estimation being ex post relative to the PAP — we acknowledge directly now.

**Not addressed:** 

### AER Ref1

> "Because the PAP is not accessible on the AEA registry, it is impossible for me to assess whether interpretation of the signals has changed after the authors obtained the results."

(Raises an implicit concern about ex-post selection of which signal to use for pinning down the structural model.)

**Addressed:** 

- PAP accessible now.

### AER Ref2

Does not mention the PAP by name, but raises a related concern about discrepancies between the AEA Registry pre-registration and the paper:

- The registry lists 150 treatment units; the paper says 144.
- The registry describes the distance variation as a *sampling* decision, not a randomized treatment arm; the paper presents it as the latter.

These discrepancies suggest either unacknowledged deviations from the registered design or unexplained discrepancies between the pre-registered and implemented protocol.

**Addressed:** 

- The 150 vs. 144 discrepancy is explained in a footnote (8 fallbacks selected; 14 clusters dropped for specific, documented field reasons).

- Treatment cell Ns by Close/Far are given for each incentive arm. 

**Not addressed:** 

- the registry framing of distance as a sampling decision vs. a randomized arm — do we address this? 

---

## Summary Table

| Source | Sample Issues | PAP / Registry Issues |
|--------|---------------|----------------------|
| QJE R1 | Complex protocol; unclear sample overlaps (take-up / baseline / observability / endline); accept-reject algorithm confusing | Structural estimation not pre-specified; framing diverges from PAP; PAP text irreconcilable with paper |
| QJE R2 | None | None |
| QJE R3 | Cluster definition unclear (158 clusters); anchor school concept unclear; centroid vs. PoT confusion; data dropped between 2.5–3.5 km | None |
| AER DLTR | None | None |
| AER Ref1 | Same cluster/PoT/centroid confusion as QJE R3; baseline N discrepancy (2,160 vs. 4,823) | PAP not accessible on AEA registry |
| AER Ref2 | 150 vs. 144 treatment units; distance described as sampling decision in registry but randomized arm in paper; "monitored at PoT" implies selection; no treatment cell Ns | Registry vs. paper discrepancies imply unacknowledged design deviation |
| AER Ref3 | Same cluster/PoT confusion; data dropped 2.5–3.5 km; 1.25 km cutoff appears arbitrary | None |
| AER Ref4 | None | None |
