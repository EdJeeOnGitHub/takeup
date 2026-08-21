#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (!length(hit)) default else substring(hit[[1L]], nchar(name) + 2L)
}
stage <- value("--stage")
if (is.null(stage) || !dir.exists(stage)) stop("--stage must name a staged candidate.")
files <- c(file.path(stage, "ECM ReStud.tex"),
           list.files(file.path(stage, "appendix", "structural-robustness", "sections"),
                      "[.]tex$", full.names = TRUE))
number_pattern <- paste0(
  "(?<![A-Za-z])(?:",
  "[0-9]+(?:[.][0-9]+)?(?:--|–|‑|-)[0-9]+(?:[.][0-9]+)?(?:\\\\%|%|\\s*percentage(?:-|‑| )points?|\\s*pp)|",
  "[0-9]+(?:[.][0-9]+)?(?:\\\\%|%|\\s*percentage(?:-|‑| )points?|\\s*pp)|",
  "(?:0|[1-9][0-9]*)[.][0-9]+",
  ")"
)
rows <- list()
for (path in files[file.exists(files)]) {
  lines <- readLines(path, warn = FALSE)
  result_language <- grepl(
    "p[[:space:]]*[=<>]|p-value|percent|percentage|[[:space:]]pp|estimate|interval|confidence|credible",
    lines, ignore.case = TRUE
  )
  relevant <- grepl("close|far|distance|signal|observab|take-up|takeup|multiplier|site", lines,
                    ignore.case = TRUE) & result_language & grepl("[0-9]", lines)
  for (line_no in which(relevant)) {
    matches <- regmatches(lines[[line_no]], gregexpr(number_pattern, lines[[line_no]], perl = TRUE))[[1L]]
    matches <- unique(matches[nzchar(matches)])
    if (!length(matches)) next
    for (old in matches) rows[[length(rows) + 1L]] <- data.frame(
      claim_id = sprintf("numeric-%04d", length(rows) + 1L),
      file = sub(paste0("^", normalizePath(stage), "/?"), "", normalizePath(path)),
      line = line_no, old_value = old, candidate_value = "",
      supporting_artifact = "", status = "TODO",
      context = trimws(lines[[line_no]]), stringsAsFactors = FALSE
    )
  }
}
audit <- if (length(rows)) do.call(rbind, rows) else data.frame()
utils::write.csv(audit, file.path(stage, "numeric-text-audit.csv"), row.names = FALSE, na = "")
header <- c("# Candidate numerical prose TODOs", "",
  "The source manuscript has not been edited. Review every item against its assigned-distance generated result.", "")
items <- if (nrow(audit)) sprintf(
  "- [ ] `%s:%d` — verify `%s` and record its candidate value/source.\n  - %s",
  audit$file, audit$line, audit$old_value, audit$context
) else "No candidate numerical claims were detected."
writeLines(c(header, items), file.path(stage, "NUMERIC-TODOS.md"))
cat(nrow(audit), "numeric prose TODOs\n")
