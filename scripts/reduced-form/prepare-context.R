#!/usr/bin/env Rscript

library(tidyverse)

script_options <- docopt::docopt(
  "Usage:
  prepare-context.R [options]

Options:
  --distance-definition=<spec>  Close/Far definition [default: assigned]
  --output-path=<path>           Output directory [default: build/assigned/context]
",
  args = if (interactive()) "" else commandArgs(trailingOnly = TRUE)
)

source("R/reduced-form/context.R")

context <- takeup_build_analysis_context(
  specification = script_options$distance_definition
)
paths <- takeup_write_analysis_context(context, script_options$output_path)
cat(paste(paths, collapse = "\n"), "\n")
