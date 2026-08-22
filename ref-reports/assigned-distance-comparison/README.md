# Assigned-versus-realized paper audit

This review report reads the active `ECM ReStud.tex`, preserves its table and
figure order, and compares regenerated realized-distance outputs on the left
with assigned-distance outputs on the right. It never promotes results into
the paper.

Run the complete local workflow with:

```sh
make distance-comparison TAKEUP_THREADS=8
```

The PDF and machine-readable artifact manifest are written to
`build/reports/assigned-vs-realized-results.pdf` and
`build/reports/assigned-vs-realized-results-manifest.csv`. Override the paper
source with `PAPER_TEX=/path/to/ECM\ ReStud.tex` when needed.

Missing Midway-only structural robustness and policy outputs appear as labeled
placeholders. Once those outputs are imported into their normal candidate build
roots, rerunning `make distance-comparison-report` fills the corresponding
panels automatically.
