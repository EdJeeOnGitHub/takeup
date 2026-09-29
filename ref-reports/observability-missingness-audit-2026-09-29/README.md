# Observability missingness audit — 29 September 2026

Initial investigation used read-only checks of saved data, generating code, and manuscript history. Subsequently applied the user-approved text edits: main-text footnote explaining 252/126, revised missingness-table sample description, and PAP count 255 to 252. The user explicitly retained the sample-construction paragraph and original closing table-note explanation. No analysis changes made. Findings below record the audit, including unresolved issues in the retained wording.

## Verified current counts

`check-counts.R` reproduces `counts.txt` from the saved clean endline and knowledge-module data.

| Current no-SMS endline sample | Count |
|---|---:|
| All clean respondents | 2,659 |
| Individual-format knowledge record (Table A) | 1,204 |
| Paired-format knowledge record (Table B) | 1,203 |
| Neither record | 252 |
| Table A respondents with positive recognition, used in analysis | 1,141 |
| Table A respondents without positive recognition | 63 |

All 252 missing no-SMS respondents have missing `survey.type`. Across the complete saved clean endline data, there are 3,678 respondents, of whom 379 have no knowledge record: 252 SMS control, 116 Social Information, 11 Reminder Only. All 379 have missing `survey.type`. These are counts in the saved clean analysis data, not a reconstruction of the manuscript's 3,729 field-completion count.

The current paper attrition generator, `scripts/appendix/paper-knowledge-table-attrition-tables.R`, asserts 1,141 observed Table A respondents and 252 missing respondents. `fit_one_sim()` samples exactly half of the 252 without replacement, producing 126 imputed missing Table A respondents per simulation (1,267 observations). It does this for 1,000 simulations. Thus 126 is not an observed assignment count. Assignment here is between individual and paired elicitation formats, not observability versus perceived observability outcomes. The code draws half globally; it does not recover each respondent's historical assignment.

## Provenance of 255

255 already appears in the first checked-in sample diagram (takeup commit `4b3c59d`, 20 February 2026), labeled "Missing knowledge module". Later diagrams and sample writeups retain it.

Crucially, `archive/code/scratch/bal-attrit-results.Rmd:84` explicitly describes 255 as missing knowledge entries among the **2,659 no-SMS endline respondents**. The same file subsequently describes the alternative sample using 252 missing respondents. This is evidence of stale/inconsistent documentation, not support for interpreting 255 as the full-SMS-inclusive count.

The exact original calculation yielding 255 has not been recovered. Current saved data give 252 under both survey-key/raw-presence and person-key/clean-summary definitions. Counting completed raw submissions linked to current no-SMS respondent IDs gives 254 missing submissions for 252 people, so duplicate submissions in this universe do not reproduce 255 either. The 2022 legacy `analysis.RData` was inspected too: its older knowledge table produces 281 missing no-SMS respondents, not 255. Do not invent an explanation for the three-person difference.

## Provenance of 2,039 and 1,784

Overleaf commit `fb8951b` (13 May 2026) introduced 2,039 into both the sample paragraph and flow diagram. Its parent described a randomly assigned subset, 1,784 module records, and 255 missing entries separately. The new version asserted 2,039 assigned and 1,784 completed plus 255 missing. Since 2,039 = 1,784 + 255, this strongly suggests an editorial sum; no independent assignment-roster calculation was located.

1,784 is reproducible: it is the number of Table A entries in `clean-endline-know-table-data.rds`. However, only 1,647 of these entries match respondents in the saved clean endline data; 137 do not. Likewise, 1,812 Table B entries comprise 1,652 matched and 160 unmatched. Therefore, 1,784 is not the count of Table A completers within the 3,678-person clean endline sample. Nor does excluding SMS alone explain its reduction to 1,141: the analysis also requires matching the analysis frame and positive recognition.

## Implications for the manuscript (not applied)

- Replace the claim that 126 observed respondents were assigned the module and missed it with the observed 252 missing no-SMS knowledge records and the half-sample reconstruction used for Table A attrition.
- In the table note, distinguish individual versus paired formats and observed versus imputed counts.
- Remove the unsupported 2,039 assignment claim; do not simply substitute 379 into prose using the different 3,729 field-completion universe.
- Reconcile the sample paragraph, flow diagram, and PAP paragraph together. Explain that 1,141 is the usable no-SMS individual-format outcome sample, not merely 1,784 minus SMS recipients.
- If the manuscript retains 3,729 and 1,784 as field/module totals, separately audit their exact universes and linkages before describing a single nesting flow.

This audit found a reporting problem; it does not change the samples or estimates in the recently rerun tables.

## Approved diagram and paragraph update

The user subsequently approved the simple diagram branch from the 2,659-person no-SMS endline sample to the 1,141-person observability analysis sample, labeled “Individual observability measure available.” Applied this to the active manuscript diagram and replaced the sample-construction paragraph with the detailed 1,204 individual / 1,203 paired / 252 missing / 1,141 usable accounting. The obsolete 2,039/1,784/255 branch was removed. The existing missingness-table closing explanation remains as previously requested. The upstream 3,729 field-completion versus 3,678 clean-record difference remains unresolved. PDF compilation succeeded and the full diagram was visually checked.

## Follow-up: endline completion reconciliation

**Current sample fully reproduced from raw records.** `reconcile-endline-completions.R` applies the production date parsing and selection rules without modifying data or rerunning outcomes. All 3,678 reconstructed respondent IDs and their selected survey keys match the saved clean dataset exactly. See `completion-flow.csv`, `contact-status-counts.csv`, `completion-checks.txt`, and the hashed record-level audit `completion-record-audit.csv`.

| Stage | Survey records | Unique people |
|---|---:|---:|
| Eligible reached/contact records after fieldwork validity, date, and GPS filters | 3,817 | 3,776 |
| Affirmative interview and consent | 3,700 | 3,678 |
| Keep first completed record per person | 3,678 | 3,678 |
| No-SMS clean sample | 2,659 | 2,659 |

The 117 records removed at the interview/consent step comprise 58 with interview=0, 21 with interview=1 and consent=0, 27 with interview=1 and consent missing, and 11 with interview missing. These are record-level categories, not mutually exclusive person-level counts. Among the 3,776 distinct contacted people, 98 have no eligible completed interview. The remaining 22 removals are repeat completed submissions. Exact status combinations are in `contact-status-counts.csv`.

The saved SMS completion counts are 163 Reminder Only + 856 Social Information = 1,019; excluding them from 3,678 yields 2,659. Thus the upstream and downstream analysis samples now reconcile exactly.

### Why the manuscript's 3,729 cannot be used as a verified starting population

The older notes report 3,830 reached (2,774 no-SMS + 1,056 SMS), minus 78 unable and 23 refusals, giving 3,729 completed (2,686 no-SMS + 1,043 SMS). Those totals are not reproduced by the production pipeline. In particular, 2,774 is the number of eligible no-SMS contact **records**, not unique people (2,745); 1,043 is the number of eligible SMS contact **records**, not completed unique SMS respondents (1,019). This exposes a unit/stage inconsistency in the old narrative, although the exact historical calculation for 1,056 and 2,686 was not found. `scripts/balance/run.R` retains the comment “Reached should be: 3830” next to the old component counts; the older sample discussion contains crossed-out/revised intermediate figures.

Do not describe the 51 difference as 51 identifiable completed respondents excluded from analysis: there is no verified 3,729-person completed roster from which to subtract them. It is the difference between an unsupported manuscript total and a fully reproduced clean total. Algebraically it combines differing reached totals, differing noncompletion totals, and duplicate handling, but that arithmetic is not a record-level exclusion explanation.

### Proposed correction, not applied

Use **3,776 distinct reached respondents in the cleaned survey records**, **98 without an eligible completed interview**, and **3,678 completed respondents**, followed by exclusion of **1,019 enrolled SMS recipients** to obtain **2,659 no-SMS respondents**. Do not retain “78 unable / 23 refused” as a breakdown of the new 98 without checking unique-person reasons; some records lack interview or consent status. State the cleaned-record scope rather than presenting 3,776 as a verified count of every field contact before validity/GPS exclusions.

Update the design data paragraph, sample-construction paragraph, PAP backup-survey account, and upstream sample-diagram boxes consistently after approval. The 4,220 field-roster count was not reaudited here. No analysis sample or estimate changes are required by this reconciliation. Manuscript remains unchanged during this investigation.

## Follow-up: refused versus unable

Checked the deployed `Endline Survey V5.xlsx` instrument and all 103 eligible records for the 98 unique noncompleters. Variables exist:

- `present`: “Is the person present?”
- `interview`: “Can the person speak to you?” (1 yes, 0 no).
- `no_interview`, with `other_no_interview`: reason unable to speak, including busy, language barrier, intoxication, underage, communication disability, mental illness, or other.
- `consent`: agreement to participate (1 yes, 0 no).
- `comments`: free-text field notes.

Unlike `present`, `interview` and `consent` are not marked required in the deployed form. The former is relevant when present=1; the latter when recruit=1. Missing fields therefore should not automatically be interpreted as refusal or inability.

`classify-noncompletion.R` gives mutually exclusive person counts, prioritizing recorded non-consent when flags overlap:

| Recorded status | People |
|---|---:|
| Any consent=0 record | 23 |
| Any interview=0 record, without consent=0 | 46 |
| Neither negative status recorded; incomplete fields | 29 |
| Total | 98 |

Three of the 23 have both non-consent and inability flags across their records. Five people have two eligible records; this is why raw counts must be deduplicated. The incomplete 29 comprise 18 with interview=1 and consent missing, 10 with consent=1 and interview missing, and one with both missing.

Free-text review shows actual reasons cross structured categories. Some interview=0/other entries explicitly describe refusal. Some consent=0 comments instead describe illness, underage status, or inability to participate. Among incomplete fields there are additional explicit refusals, scheduling issues, absence, and accidental/premature form finalizations. Six of the 29 have no comment. Thus “23 recorded non-consent” is reproducible; “exactly 23 refusals and 75 unable” is not a validated substantive classification. A strict two-way refused/unable split would require adjudication and still needs an unknown/administrative category.

For concise manuscript reporting, use “98 did not complete” or, if a breakdown is wanted, “23 recorded non-consent; 75 other or incomplete records.” Do not restore the old “78 unable” count. No manuscript or analysis changes made in this follow-up.

## Approved final count update

Applied the user's “Did not complete” choice to the active manuscript: 3,776 reached, 98 did not complete, 3,678 completed, with 1,019 SMS recipients excluded to yield 2,659. The text and figure notes specify distinct respondents in the cleaned survey records. Updated the design paragraph, sample-construction paragraph, PAP account, and diagram. The source diagram in `presentations/sample-flow-diagram/` and frozen replication copy now match the active manuscript diagram; the artifact-contract hash was refreshed. Updated the source sample writeup as well. Older historical notes, archived original diagrams, and review comparisons are retained as historical records. No data or results tables changed.
