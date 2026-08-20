#!/usr/bin/env Rscript

tex_path <- Sys.getenv("TAKEUP_PAPER_TEX", "ref-reports/ECM ReStud.tex")
if (!file.exists(tex_path)) stop("Paper source not found: ", tex_path)
lines <- readLines(tex_path, warn = FALSE)
strip_comment <- function(x) sub("(?<!\\\\)%.*$", "", x, perl = TRUE)
active <- vapply(lines, strip_comment, character(1))

extract_paths <- function(pattern, kind) {
  hits <- regmatches(active, gregexpr(pattern, active, perl = TRUE))
  values <- unlist(hits, use.names = FALSE)
  if (!length(values)) return(data.frame())
  path <- sub("^[^{]+\\{", "", values)
  path <- sub("\\}$", "", path)
  data.frame(kind = kind, manuscript_path = path, stringsAsFactors = FALSE)
}

registry <- rbind(
  extract_paths("\\\\input\\{[^}]+\\}", "input"),
  extract_paths("\\\\includegraphics(?:\\[[^]]*\\])?\\{[^}]+\\}", "figure")
)
registry <- unique(registry)
registry$expected_extension <- ifelse(
  grepl("\\.[A-Za-z0-9]+$", registry$manuscript_path),
  "", ifelse(registry$kind == "input", ".tex", ".pdf")
)
registry$artifact <- paste0(registry$manuscript_path,
                            registry$expected_extension)

roots <- c(
  "presentations", "appendix/structural-robustness", "optim",
  "temp-data/distributional-robustness", "temp-data/bimodal", "temp-data", "."
)
resolve <- function(path) {
  aliases <- c(
    "praise-stigma-by-externality-knowledge.tex" =
      "presentations/rf-tables/main-specs/new-praise-stigma-by-externality-knowledge.tex"
  )
  alias <- unname(aliases[basename(path)])
  if (length(alias) && !is.na(alias) && file.exists(alias)) return(alias)
  candidates <- unique(c(path, basename(path), file.path(roots, path)))
  direct <- candidates[file.exists(candidates)]
  if (length(direct)) return(direct[[1L]])
  basename_matches <- list.files(
    ".", recursive = TRUE, full.names = TRUE,
    pattern = paste0("^", gsub("([.(){}+*?^$|\\[\\]\\\\])", "\\\\\\1",
                                basename(path)), "$")
  )
  basename_matches <- basename_matches[!grepl("^(./)?build/", basename_matches)]
  if (length(basename_matches)) basename_matches[[1L]] else NA_character_
}
registry$source_path <- vapply(registry$artifact, resolve, character(1))
registry$status <- ifelse(is.na(registry$source_path), "missing", "resolved")
registry$generator_group <- ifelse(
  grepl("rf-|rf_|reduced|belief|takeup", registry$artifact, ignore.case = TRUE),
  "reduced-form",
  ifelse(grepl("balance|attrition|implementation|travel", registry$artifact,
               ignore.case = TRUE), "balance", "other")
)

dir.create("build/manifest", recursive = TRUE, showWarnings = FALSE)
path <- "build/manifest/paper-artifacts.csv"
write.csv(registry, path, row.names = FALSE)
cat("Wrote", nrow(registry), "active manuscript dependencies to", path, "\n")
cat("Unresolved:", sum(registry$status == "missing"), "\n")
