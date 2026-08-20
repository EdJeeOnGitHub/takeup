# Reusable R code

Files below `R/` define reusable functions and constants. Sourcing them must
not launch an analysis merely because the file was loaded. Executable tasks,
argument parsing, data loading with side effects, and output writing belong in
`scripts/`.

The reduced-form and balance workflows share an explicit cached context built
by `R/reduced-form/context.R`. The module defines readers and writers, but only
an entrypoint or declared pipeline target may call the writing function.
