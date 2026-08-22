#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1L]])
}

paper_tex <- value(
  "--paper-tex",
  "/home/agent/projects/overleaf/overleaf-takeup/ECM ReStud.tex"
)
output_root <- value("--output-root", "build/reports")
input <- "ref-reports/assigned-distance-comparison/assigned-vs-realized-results.Rmd"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  input,
  output_file = "assigned-vs-realized-results.pdf",
  output_dir = normalizePath(output_root),
  knit_root_dir = normalizePath("."),
  params = list(paper_tex = paper_tex, output_root = output_root),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
