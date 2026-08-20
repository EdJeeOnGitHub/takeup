#!/usr/bin/env Rscript

modules <- c(
  "R/data/survey-cleaning.R",
  "R/reduced-form/context.R"
)

old_options <- options()
old_env <- Sys.getenv()
old_wd <- getwd()
old_files <- list.files(".", all.files = TRUE, recursive = TRUE)

module_env <- new.env(parent = globalenv())
for (module in modules) sys.source(module, envir = module_env)

stopifnot(
  identical(getwd(), old_wd),
  identical(options(), old_options),
  identical(Sys.getenv(), old_env),
  identical(list.files(".", all.files = TRUE, recursive = TRUE), old_files)
)

required_functions <- c(
  "prepare_baseline_data", "prepare_endline_data",
  "takeup_build_analysis_context", "takeup_read_analysis_context",
  "takeup_write_analysis_context"
)
stopifnot(all(vapply(required_functions, exists, logical(1),
                     envir = module_env, inherits = FALSE)))

cat("Reusable analysis-context modules source without side effects.\n")
