# Manuscript promotion — 29 September 2026

Promoted after the user authorized bringing all rerun outputs into the manuscript and requested syncing tables followed by checking text references.

- 13 affected tables: audited corrected data cells, preserving manuscript labels and layout.
- Two numerically unchanged no-controls robustness tables: corrected significance stars; observability also uses proper LaTeX less-than signs.
- Three prose files: introduction (`ECM ReStud.tex`), `experimental-design.tex`, and `online-appendix.tex`. Reconciled estimates, p-values, Lee bounds, attrition tests, and baseline-imbalance robustness statistics. Filled all three empty p-values.
- Existing interpretation retained except the module-attrition wording now distinguishes significance at 5% from the pooled Bracelet difference at 10%, and clarifies joint tests across distance cells.

`before/` preserves original files, `staged/` contains promoted files, per-file `.diff` files record exact edits, and `manifest.json` records before/after SHA256 hashes. `prepare.py` records preparation of the initial 16 files; the two star-only tables were appended afterward.

Validation: all 18 promoted file hashes match the live manuscript; all data cells in the 15 tables match the audited corrected outputs; no empty p-values remain in the updated prose. No new estimation was required. The PDF rebuild succeeds. Existing warnings remain for the undefined PAP reference `tab:endline-obs-attrit-backup` and four oversized floats; these remain separate manuscript TODOs.
