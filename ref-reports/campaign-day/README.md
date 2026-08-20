# Campaign-day original-assignment analysis

This directory is an isolated, local appendix analysis. It does not alter or
enter the main paper.

Run the estimator from the repository root in the pinned R environment:

```sh
Rscript scripts/appendix/campaign-day-original-assignment.R \
  --bootstrap-draws=1000 \
  --seed=20260807
```

Then compile the note from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error campaign-day-note.tex
```

The estimator uses directly recorded take-up and campaign day, matching the
original-assignment endpoint. It deliberately excludes IV and conditional
hazard specifications. Generated data include a day-12 reproduction audit,
complete tidy summaries, and the paired cluster-weight draws.
