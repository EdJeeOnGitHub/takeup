strip_latex_comment <- function(line) {
  sub("(?<!\\\\)%.*$", "", line, perl = TRUE)
}

extract_braced_value <- function(text, command) {
  marker <- regexpr(paste0("\\\\", command, "(?:\\[[^]]*\\])?\\{"),
                    text, perl = TRUE)
  if (marker[[1L]] < 0L) return(NA_character_)
  start <- marker[[1L]] + attr(marker, "match.length")
  chars <- strsplit(substr(text, start, nchar(text)), "", fixed = TRUE)[[1L]]
  depth <- 1L
  escaped <- FALSE
  value <- character()
  for (char in chars) {
    if (identical(char, "\\") && !escaped) {
      escaped <- TRUE
      value <- c(value, char)
      next
    }
    if (!escaped && identical(char, "{")) depth <- depth + 1L
    if (!escaped && identical(char, "}")) {
      depth <- depth - 1L
      if (depth == 0L) break
    }
    value <- c(value, char)
    escaped <- FALSE
  }
  trimws(paste0(value, collapse = ""))
}

nearest_float_caption <- function(lines, line_number) {
  before <- lines[seq_len(line_number)]
  begins <- grep("\\\\begin\\{(table|figure)\\}", before, perl = TRUE)
  ends <- grep("\\\\end\\{(table|figure)\\}", before, perl = TRUE)
  begin <- if (length(begins)) tail(begins, 1L) else NA_integer_
  end <- if (length(ends)) tail(ends, 1L) else NA_integer_
  if (is.na(begin) || (!is.na(end) && end > begin)) return(NA_character_)
  after <- grep("\\\\end\\{(table|figure)\\}",
                lines[line_number:length(lines)], perl = TRUE)
  float_end <- if (length(after)) line_number + after[[1L]] - 1L else line_number
  extract_braced_value(paste(lines[begin:float_end], collapse = "\n"), "caption")
}

nearest_float_label <- function(lines, line_number) {
  before <- lines[seq_len(line_number)]
  begins <- grep("\\\\begin\\{(table|figure)\\}", before, perl = TRUE)
  ends <- grep("\\\\end\\{(table|figure)\\}", before, perl = TRUE)
  begin <- if (length(begins)) tail(begins, 1L) else NA_integer_
  end <- if (length(ends)) tail(ends, 1L) else NA_integer_
  if (is.na(begin) || (!is.na(end) && end > begin)) return(NA_character_)
  after <- grep("\\\\end\\{(table|figure)\\}",
                lines[line_number:length(lines)], perl = TRUE)
  float_end <- if (length(after)) line_number + after[[1L]] - 1L else line_number
  extract_braced_value(paste(lines[begin:float_end], collapse = "\n"), "label")
}

nearest_section <- function(lines, line_number) {
  before <- lines[seq_len(line_number)]
  hits <- grep("\\\\(section|subsection|subsubsection)\\*?\\{", before,
               perl = TRUE)
  if (!length(hits)) return("Manuscript results")
  text <- before[tail(hits, 1L)]
  command <- sub(".*\\\\(section|subsection|subsubsection).*", "\\1", text)
  value <- extract_braced_value(text, command)
  ifelse(is.na(value) || !nzchar(value), "Manuscript results", value)
}

parse_paper_artifacts <- function(paper_tex) {
  raw <- readLines(paper_tex, warn = FALSE)
  lines <- vapply(raw, strip_latex_comment, character(1))
  pattern <- paste0(
    "\\\\input\\{[^}]+\\}|",
    "\\\\includegraphics(?:\\[[^]]*\\])?\\{[^}]+\\}"
  )
  rows <- list()
  order <- 0L
  for (line_number in seq_along(lines)) {
    matches <- regmatches(lines[[line_number]],
                          gregexpr(pattern, lines[[line_number]], perl = TRUE))[[1L]]
    if (!length(matches) || identical(matches, "")) next
    for (match in matches) {
      order <- order + 1L
      kind <- if (startsWith(match, "\\input")) "table" else "figure"
      path <- sub("^[^{]+\\{", "", match)
      path <- sub("\\}$", "", path)
      extension <- if (grepl("\\.[A-Za-z0-9]+$", path)) "" else
        if (kind == "table") ".tex" else ".pdf"
      caption <- nearest_float_caption(lines, line_number)
      label <- nearest_float_label(lines, line_number)
      included_source <- file.path(dirname(paper_tex), paste0(path, extension))
      if ((is.na(caption) || is.na(label)) && file.exists(included_source) &&
          grepl("[.]tex$", included_source)) {
        included_text <- paste(readLines(included_source, warn = FALSE),
                               collapse = "\n")
        if (is.na(caption)) caption <- extract_braced_value(included_text, "caption")
        if (is.na(label)) label <- extract_braced_value(included_text, "label")
      }
      rows[[length(rows) + 1L]] <- data.frame(
        paper_order = order,
        paper_line = line_number,
        kind = kind,
        artifact = paste0(path, extension),
        section = nearest_section(lines, line_number),
        caption = caption,
        label = label,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

read_overleaf_labels <- function(aux_path) {
  if (!file.exists(aux_path)) {
    return(data.frame(label = character(), overleaf_number = character()))
  }
  lines <- readLines(aux_path, warn = FALSE)
  pattern <- "^\\\\newlabel\\{([^}]+)\\}\\{\\{([^}]*)\\}"
  matches <- regexec(pattern, lines, perl = TRUE)
  values <- regmatches(lines, matches)
  values <- values[lengths(values) == 3L]
  if (!length(values)) {
    return(data.frame(label = character(), overleaf_number = character()))
  }
  data.frame(
    label = vapply(values, `[[`, character(1L), 2L),
    overleaf_number = vapply(values, `[[`, character(1L), 3L),
    stringsAsFactors = FALSE
  )
}

infer_paper_numbers <- function(paper_tex) {
  lines <- vapply(readLines(paper_tex, warn = FALSE), strip_latex_comment,
                  character(1L))
  appendix <- FALSE
  subsection <- 0L
  counters <- c(table = 0L, figure = 0L)
  rows <- list()
  index <- 1L
  while (index <= length(lines)) {
    line <- lines[[index]]
    if (grepl("^\\s*\\\\startonlineappendix\\s*$", line, perl = TRUE)) {
      appendix <- TRUE
      subsection <- 0L
      counters[] <- 0L
    }
    if (appendix && grepl("\\\\subsection\\{", line, perl = TRUE)) {
      subsection <- subsection + 1L
      counters[] <- 0L
    }
    begin <- regexec("\\\\begin\\{(table|figure)\\}", line, perl = TRUE)
    value <- regmatches(line, begin)[[1L]]
    if (length(value) == 2L) {
      kind <- value[[2L]]
      counters[[kind]] <- counters[[kind]] + 1L
      end_pattern <- paste0("\\\\end\\{", kind, "\\}")
      end <- index
      while (end <= length(lines) &&
             !grepl(end_pattern, lines[[end]], perl = TRUE)) end <- end + 1L
      block <- paste(lines[index:min(end, length(lines))], collapse = "\n")
      label <- extract_braced_value(block, "label")
      if (!is.na(label) && nzchar(label)) {
        number <- if (appendix) {
          paste0(intToUtf8(64L + subsection), counters[[kind]])
        } else {
          as.character(counters[[kind]])
        }
        rows[[length(rows) + 1L]] <- data.frame(
          label = label, overleaf_number = number, stringsAsFactors = FALSE
        )
      }
      index <- end
    }
    index <- index + 1L
  }
  if (!length(rows)) {
    return(data.frame(label = character(), overleaf_number = character()))
  }
  do.call(rbind, rows)
}

read_artifact_contract <- function(path = "replication/paper-artifact-contract.csv") {
  contract <- utils::read.csv(path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  contract[, c("artifact", "workflow_owner", "default_contract")]
}

match_contract <- function(artifact, contract) {
  exact <- which(contract$artifact == artifact)
  if (length(exact) == 1L) return(contract[exact, , drop = FALSE])
  base <- which(basename(contract$artifact) == basename(artifact))
  if (length(base) == 1L) return(contract[base, , drop = FALSE])
  data.frame(
    artifact = artifact,
    workflow_owner = "unclassified",
    default_contract = "unknown",
    stringsAsFactors = FALSE
  )
}

is_distance_sensitive <- function(owner, contract) {
  identical(contract, "generated") || owner %in% c(
    "reduced_form", "balance", "structural_postprocess", "optimal_policy"
  )
}

resolve_artifact <- function(artifact, roots) {
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(NA_character_)
  direct <- file.path(roots, artifact)
  direct <- direct[file.exists(direct)]
  if (length(direct)) return(normalizePath(direct[[1L]]))
  suffix_hits <- unlist(lapply(roots, function(root) {
    files <- list.files(root, recursive = TRUE, full.names = TRUE)
    files[endsWith(files, artifact)]
  }), use.names = FALSE)
  if (length(suffix_hits)) return(normalizePath(suffix_hits[[1L]]))
  basename_hits <- unlist(lapply(roots, function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE,
               pattern = paste0("^", gsub("([.(){}+*?^$|\\[\\]\\\\])",
                                             "\\\\\\1", basename(artifact)), "$"))
  }), use.names = FALSE)
  if (length(basename_hits)) return(normalizePath(basename_hits[[1L]]))
  NA_character_
}

sha256_file <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

table_fragment_columns <- function(path) {
  if (is.na(path) || !file.exists(path) || !grepl("[.]tex$", path)) return(NA_integer_)
  lines <- readLines(path, warn = FALSE)
  hit <- grep("\\\\begin\\{tabular\\}(?:\\[[^]]*\\])?\\{[^}]+\\}",
              lines, value = TRUE, perl = TRUE)
  if (!length(hit)) return(NA_integer_)
  specification <- sub(".*\\{([^{}]+)\\}.*", "\\1", hit[[1L]])
  lengths(regmatches(specification, gregexpr("[lcrpmbX]", specification,
                                              perl = TRUE)))
}

table_fragment_rows <- function(path) {
  if (is.na(path) || !file.exists(path) || !grepl("[.]tex$", path)) return(NA_integer_)
  lines <- readLines(path, warn = FALSE)
  sum(grepl("\\\\\\\\$", trimws(lines)))
}

build_comparison_manifest <- function(paper_tex, realized_roots,
                                      assigned_roots, contract_path) {
  artifacts <- parse_paper_artifacts(paper_tex)
  contract <- read_artifact_contract(contract_path)
  metadata <- lapply(artifacts$artifact, match_contract, contract = contract)
  artifacts$workflow_owner <- vapply(metadata, `[[`, character(1L),
                                     "workflow_owner")
  artifacts$default_contract <- vapply(metadata, `[[`, character(1L),
                                       "default_contract")
  artifacts$distance_sensitive <- mapply(
    is_distance_sensitive,
    artifacts$workflow_owner,
    artifacts$default_contract
  )
  artifacts$realized_path <- vapply(
    artifacts$artifact, resolve_artifact, character(1L), roots = realized_roots
  )
  artifacts$assigned_path <- vapply(
    artifacts$artifact, resolve_artifact, character(1L), roots = assigned_roots
  )
  artifacts$realized_sha256 <- vapply(artifacts$realized_path, sha256_file,
                                      character(1L))
  artifacts$assigned_sha256 <- vapply(artifacts$assigned_path, sha256_file,
                                      character(1L))
  artifacts$status <- ifelse(
    is.na(artifacts$realized_path) & is.na(artifacts$assigned_path),
    "both missing",
    ifelse(is.na(artifacts$realized_path), "realized missing",
      ifelse(is.na(artifacts$assigned_path), "assigned missing",
        ifelse(artifacts$realized_sha256 == artifacts$assigned_sha256,
               "identical", "changed")))
  )
  realized_columns <- vapply(artifacts$realized_path, table_fragment_columns,
                             integer(1L))
  assigned_columns <- vapply(artifacts$assigned_path, table_fragment_columns,
                             integer(1L))
  artifacts$table_columns <- pmax(realized_columns, assigned_columns,
                                  na.rm = TRUE)
  artifacts$table_columns[is.infinite(artifacts$table_columns)] <- NA_integer_
  artifacts$wide_layout <- artifacts$kind == "table" &
    !is.na(artifacts$table_columns) & artifacts$table_columns > 8L
  realized_rows <- vapply(artifacts$realized_path, table_fragment_rows,
                          integer(1L))
  assigned_rows <- vapply(artifacts$assigned_path, table_fragment_rows,
                          integer(1L))
  artifacts$table_rows <- pmax(realized_rows, assigned_rows, na.rm = TRUE)
  artifacts$table_rows[is.infinite(artifacts$table_rows)] <- NA_integer_
  artifacts$two_page_layout <- artifacts$wide_layout &
    !is.na(artifacts$table_rows) & artifacts$table_rows > 16L
  artifacts$detail_page <- artifacts$distance_sensitive &
    artifacts$status != "identical"
  artifacts
}

latex_path <- function(path) {
  paste0("\\detokenize{", normalizePath(path, mustWork = FALSE), "}")
}

emit_placeholder <- function(label, artifact, owner) {
  cat(
    "\\fcolorbox{gray}{gray!8}{%\n",
    "\\begin{minipage}[c][0.48\\textheight][c]{0.92\\linewidth}\n",
    "\\centering\\large ", label, " unavailable\\par\\vspace{1em}\n",
    "\\small\\texttt{\\detokenize{", artifact, "}}\\par\\vspace{0.6em}\n",
    "\\small Expected from workflow: ",
    paste(strsplit(owner, "_", fixed = TRUE)[[1L]], collapse = "\\_"), "\n",
    "\\end{minipage}}\n",
    sep = ""
  )
}

emit_table <- function(path) {
  text <- readLines(path, warn = FALSE)
  text <- text[!grepl("^\\s*\\\\(begin|end)\\{table\\}", text)]
  cat(text, sep = "\n")
  cat("\n")
}

emit_figure <- function(path) {
  cat("\\includegraphics[width=\\linewidth,height=0.70\\textheight,keepaspectratio]{",
      latex_path(path), "}\n", sep = "")
}

emit_panel <- function(label, path, kind, artifact, owner,
                       width = "0.49\\linewidth") {
  cat("\\begin{minipage}[t]{", width, "}\n\\centering\n", sep = "")
  cat("\\textbf{", label, "}\\par\\vspace{0.35em}\n", sep = "")
  if (is.na(path)) {
    emit_placeholder(label, artifact, owner)
  } else if (kind == "table") {
    cat("\\begingroup\\scriptsize\n")
    emit_table(path)
    cat("\\endgroup\n")
  } else {
    emit_figure(path)
  }
  cat("\\vspace{0.4em}\\par{\\tiny\\texttt{\\detokenize{",
      ifelse(is.na(path), artifact, path), "}}}\n", sep = "")
  cat("\\end{minipage}\n")
}

emit_comparison_page <- function(row) {
  caption <- row$caption
  if (is.na(caption) || !nzchar(caption)) {
    caption <- paste0("\\texttt{\\detokenize{", basename(row$artifact), "}}")
  }
  numbered_caption <- if (!is.na(row$overleaf_number) &&
                          nzchar(row$overleaf_number)) {
    paste0(ifelse(row$kind == "table", "Table ", "Figure "),
           row$overleaf_number, ": ", caption)
  } else {
    caption
  }
  cat("\\section*{", row$section, "}\n", sep = "")
  cat("\\begin{center}\\large\\textbf{", numbered_caption, "}\\end{center}\n",
      sep = "")
  if (isTRUE(row$two_page_layout)) {
    emit_panel("Realized distance", row$realized_path,
               row$kind, row$artifact, row$workflow_owner,
               width = "\\linewidth")
    cat("\\par\\vspace{0.5em}{\\footnotesize Assigned-distance comparison continues on the next page.}\\par\n")
    cat("\\clearpage\n")
    cat("\\begin{center}\\large\\textbf{", numbered_caption,
        " (continued)}\\end{center}\n", sep = "")
    emit_panel("Assigned distance", row$assigned_path,
               row$kind, row$artifact, row$workflow_owner,
               width = "\\linewidth")
  } else if (isTRUE(row$wide_layout)) {
    emit_panel("Realized distance", row$realized_path,
               row$kind, row$artifact, row$workflow_owner,
               width = "\\linewidth")
    cat("\\par\\vspace{0.8em}\n")
    emit_panel("Assigned distance", row$assigned_path,
               row$kind, row$artifact, row$workflow_owner,
               width = "\\linewidth")
  } else {
    emit_panel("Realized distance", row$realized_path,
               row$kind, row$artifact, row$workflow_owner)
    cat("\\hfill\n")
    emit_panel("Assigned distance", row$assigned_path,
               row$kind, row$artifact, row$workflow_owner)
  }
  cat("\\par\\vspace{0.5em}{\\footnotesize Status: ", row$status,
      "; manuscript line ", row$paper_line, ".}\\par\n", sep = "")
}
