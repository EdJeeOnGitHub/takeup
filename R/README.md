# Reusable R code

Files below `R/` define reusable functions and constants. Sourcing them must
not launch an analysis merely because the file was loaded. Executable tasks,
argument parsing, data loading with side effects, and output writing belong in
`scripts/`.
