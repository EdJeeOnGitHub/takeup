# Sample Flow Diagram — Plan

## Goal

A single-page TikZ figure for the working paper. Clean, consistent taxonomy.

---

## What economists/trialists usually do

CONSORT-style flow diagrams (JAMA, AER, JPE, etc.) follow a standard grammar:

| Stage | Label convention | Visual convention |
|---|---|---|
| Population / frame | "Assessed for eligibility" | Plain box |
| Exclusions | "Excluded (reason: N)" | Indented or side box |
| Randomised | "Randomised (N)" | Box + horizontal split into arms |
| Allocated | "Allocated to [arm] (N)" | Column header box |
| Lost to follow-up | "Lost to follow-up (N)" | Side branch |
| Analysed | "Analysed (N)" | Shaded box at bottom |

Key principles:
- One column per arm, parallel vertical flow.
- Side branches for exclusions/attrition only — not for design notes.
- Footnotes carry implementation detail; the diagram shows counts only.
- Two visual dimensions at most (fill color + border style).

---

## Proposed taxonomy for this study

### Columns: data sources, not treatment arms

| Column | Label | What it tracks |
|---|---|---|
| A | **Baseline** | Adults surveyed before treatment |
| B | **Monitored sample** | Adults observed for deworming take-up at PoT |
| C | **Endline** | Adults surveyed after treatment |

This organizes the diagram around *what we measured*, not *who got which SMS*.
Treatment arm (SMS treatment / SMS control) is an attribute of individuals
within columns B and C, not a column in itself. Carry treatment-arm breakdowns
as sub-counts within boxes or in footnotes.

### Stage labels per column

**A — Baseline**
1. Selected from census roster
2. Targeted for interview
3. Reached
4. Completed baseline *(analysis sample)*
5. Not completed *(side branch)*

**B — Monitored sample**
1. Three contributing sub-lists, each drawn from census:
   - Baseline roster (4,185)
   - SMS treatment enrolled (3,022)
   - SMS control roster (7,155: 5,803 phone + 1,352 non-phone)
2. Total monitored, unique individuals (deduplicated)
3. Main take-up analysis sample *(excl. SMS treatment)*

**C — Endline**
1. Sampling frame (SMS treatment enrolled + SMS control roster)
2. Selected for endline
3. Reached
4. Completed endline *(analysis sample)*
5. Missing knowledge module *(side branch, deviation)*

### Visual styles (two dimensions only)

| Style | Meaning | Implementation |
|---|---|---|
| Plain box | Roster / field count | `draw, rounded corners` |
| Green fill | Completed survey (analysis sample) | `fill=green!15` |
| Blue fill | Monitored at PoT (take-up observed) | `fill=blue!10` |
| Orange fill | Known deviation / data issue | `fill=orange!20` |
| Dashed border | Derived / explanatory note | `dashed` |

Drop the PAP-target and code-note boxes entirely; move to footnotes.

---

## Proposed simplified counts

(From the Mermaid prototype; verify against `prepare_analysis_data.R` output.)

**Shared top: Census frame**
- 144 clusters, 39,301 adults (18+)

**Column A — Baseline**
- Selected: 4,185 adults (2,414 HH)
- Targeted for interview: 2,160
- Reached: 2,062
- Completed: 1,995 *(green)*
- Not completed: 67 (refused 29; unable 38)

**Column B — Monitored sample**
- Baseline roster: 4,185
- SMS treatment (enrolled): 3,022 (SI: 2,635; reminder: 387)
- SMS control: 7,155 (phone-owners 5,803; non-phone 1,352)
  - Note: 989 non-phone individuals dropped due to SurveyCTO error *(orange)*
- Total unique, deduplicated: 12,827 *(blue)*
- Main take-up analysis (excl. SMS treatment): 9,805 *(blue, dashed border)*

**Column C — Endline**
- Sampling frame: 11,166 (SMS treatment 3,022 + SMS control 8,144)
- Selected: 4,220 (3,586 main + 634 backups)
- Reached: 3,830 (SMS treatment 1,056; SMS control 2,774)
- Completed: 3,729 *(green)* (SMS treatment 1,043; SMS control 2,686)
- Knowledge module missing: 255 *(orange)*

---

## Decisions (resolved)

1. **Column order**: Baseline → Monitored sample → Endline (A–B–C).
2. **Monitored sample feeders**: three separate feeder boxes (Baseline roster,
   SMS treatment, SMS control), each flowing down into the merged total.
3. **989-drop deviation**: shown as an orange side box off the SMS control feeder.
4. **Sub-count breakdowns**: totals only in each box; no arm-level breakdowns.

---

## Layout sketch

```
  ┌─────────────────────────────────────────────────────────────────────┐
  │         Census frame: 144 clusters, 39,301 adults (18+)             │
  └────┬──────────────────────────────────────────────────────────┬─────┘
       │                                                          │
  [A: Baseline]            [B: Monitored sample]          [C: Endline]
       │                                                          │
  Selected            ┌──────────┬────────────┐          Sampling frame
   4,185              │          │            │              11,166
       │         Baseline   SMS Tx       SMS Ctrl               │
  Targeted         4,185     3,022         7,155           Selected
   2,160              │          │            │  \           4,220
       │              │          │            │  989 🟠        │
  Reached             └──────────┴────────────┘           Reached
   2,062                         │                          3,830
       │                  Total unique                         │
  Completed               12,827 🔵                      Completed
   1,995 🟢                      │                        3,729 🟢
       │                 Take-up analysis                       │
  Not completed          9,805 (dashed) 🔵              Missing K-table
      67                                                   255 🟠
```

---

## Implementation plan for `sample-consort.tex`

1. Use `tikzpicture` with `matrix of nodes` or manual `\node` placement.
   - Manual placement (like the LaTeX prototype) gives more control.
   - Use `node distance` relative coordinates.

2. Define colors in preamble:
   ```latex
   \colorlet{analysisgreen}{green!15}
   \colorlet{monblue}{blue!10}
   \colorlet{issueoran}{orange!20}
   ```

3. Three-column layout with shared top (census) and shared bottom (dedup box).

4. Arrow styles:
   - `->` solid for main flow
   - `->, dashed` for derived/note connections
   - Side branches use `|-` or `to[out=..., in=...]` bends

5. Group boxes (subgraph equivalents): `\node[fit=..., draw, rounded corners]`
   with a label above.

6. Legend: small inset `tikzpicture` or a `tabular` below the figure.

7. Target: fits on one A4/letter page. Use `\resizebox{\textwidth}{!}{...}`.

---

## Decisions needed from you

- [ ] Split arm C into phone/non-phone sub-columns, or keep one column?
- [ ] Endline rows appended to columns B+C, or a separate merged section?
- [ ] Include the dedup/aggregate box in the figure, or just report in text?
- [ ] Any arms or stages to drop entirely (e.g., drop baseline monitoring)?
