# Final PDF inspection — 29 September 2026

Inspected the current 122-page manuscript PDF: visual overview of every page via contact sheets, then full-size checks of pages 64, 88, 96, 115, 118 and 122. Reviewed the build log and extracted text. This is a layout/reference inspection, not a fresh audit of all estimates. No manuscript source changes made during this inspection. The alternative-model review was marked resolved at the user's direction.

## Findings

- No undefined-reference or citation warnings in the current build log. Main-text layout, main regression tables, and revised sample diagram (page 88) have no obvious clipping or overlapping elements.
- **Page 64, Table B15:** three-panel knowledge/accuracy table is too tall; notes collide with the page number. Log reports 42.52 pt excess height.
- **Page 115, Table L1:** structural robustness table notes collide with the page number. Log reports 29.65 pt excess height.
- **Page 118, Table M1:** distance-cap table notes cross the page-number area and approach the page bottom. Log reports 56.90 pt excess height.
- **Page 117:** appendix M heading and subsection heading occupy an otherwise empty page before Table M1. Keep the heading with its table when adjusting that float.
- **Page 122:** the final inline formula protrudes into the right margin (40.56 pt overfull box). A displayed equation or explicit break would fix this without changing the mathematics.
- **Page 96:** “Additions beyond the PAP” is stranded at the bottom of the continued PAP table; its entries start on page 97. Keep the heading with the next row. The PAP table also has 32 pt width warnings, though its content remains inside the physical page.

There are **three**, not four, remaining oversized-float warnings in this build. Other small overfull boxes exist; no obvious physical-page clipping was apparent in the overview. Many appendix floats have generous whitespace, but that alone is not an error.

## Recommended layout work

Reduce height/spacing of B15, L1, and M1 while keeping notes readable and within the text area; keep the appendix M heading with M1; move the PAP group heading to its first row's page; display/break the final formula. Rebuild and visually recheck those pages before calling layout finalized.

## Missingness audit still pending

The approved observability branch is reconciled: 2,659 no-SMS respondents = 1,204 individual + 1,203 paired + 252 missing; 1,141 of the individual-format respondents recognize at least one peer. The attrition analysis imputes 126 missing individual-format respondents per simulation.

What remains is tracing the manuscript's 3,729 completed surveys to the 3,678 clean endline respondents (a difference of 51). This requires establishing the source of the field-completion count and the record-level exclusions/deduplication/linkage bridge; the difference must not be attributed to any one cleaning step without evidence. Then reconcile the upstream diagram and the corresponding design, sample-construction and PAP passages. No new estimation has yet been shown necessary.

The historical source calculation for 255 remains unknown, but that count has been removed from the approved appendix paragraph and active sample diagram, so recovering it is a historical audit question rather than a prerequisite for those edits. The original table-note format terminology remains as explicitly requested by the user.

PDF SHA256: `dc0e071c704382a95becc2f4a5096453510f3549ec6213863487cca014fcbf77`
