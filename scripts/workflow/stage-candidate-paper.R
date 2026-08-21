#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (!length(hit)) default else substring(hit[[1L]], nchar(name) + 2L)
}
specification <- value("--spec", Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned"))
source("R/distance/spec.R")
specification <- takeup_distance_spec(specification)
stage <- value("--stage", file.path("build", "paper-candidate", specification))
strict <- "--strict" %in% args

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  sub("[[:space:]].*$", "", system2("sha256sum", path, stdout = TRUE)[[1L]])
}
strip_comment <- function(x) sub("(?<!\\\\)%.*$", "", x, perl = TRUE)
escape_regex <- function(x) gsub("([.(){}+*?^$|\\[\\]\\\\])", "\\\\\\1", x)
copy_one <- function(source, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  isTRUE(file.copy(source, destination, overwrite = TRUE))
}

status <- system2("Rscript", c("--no-save", "--no-restore",
  "scripts/workflow/build-paper-artifact-registry.R"))
if (!identical(status, 0L)) stop("Could not build main-paper registry.")
main <- read.csv("build/manifest/paper-artifacts.csv", stringsAsFactors = FALSE)
contract <- read.csv("replication/paper-artifact-contract.csv", stringsAsFactors = FALSE,
                     na.strings = character())
main <- merge(main, contract, by = "artifact", all.x = TRUE, sort = FALSE,
              suffixes = c("", "_contract"))
main$scope <- "manuscript"
main$current_path <- ifelse(
  main$default_contract == "frozen", main$deposit_path, main$source_path
)

appendix_root <- "appendix/structural-robustness"
appendix_tex <- c(file.path(appendix_root, "structural-robustness-appendix.tex"),
                  list.files(file.path(appendix_root, "sections"), "[.]tex$",
                             full.names = TRUE))
extract_appendix <- function(path) {
  text <- paste(strip_comment(readLines(path, warn = FALSE)), collapse = "\n")
  text <- gsub("\\\\StructuralRobustnessPath", appendix_root, text)
  patterns <- c("\\\\input\\{([^}]+)\\}",
                "\\\\includegraphics(?:\\[[^]]*\\])?\\{([^}]+)\\}")
  values <- unlist(lapply(patterns, function(pattern) {
    hits <- gregexpr(pattern, text, perl = TRUE)
    matches <- regmatches(text, hits)[[1L]]
    if (!length(matches) || identical(matches, character())) return(character())
    sub(pattern, "\\1", matches, perl = TRUE)
  }), use.names = FALSE)
  values <- values[grepl(paste0("^", escape_regex(appendix_root), "/(tables|figures)/"), values)]
  values <- ifelse(grepl("[.][A-Za-z0-9]+$", values), values,
                   paste0(values, ifelse(grepl("/tables/", values), ".tex", ".pdf")))
  unique(values)
}
appendix_artifacts <- unique(unlist(lapply(appendix_tex, extract_appendix), use.names = FALSE))
appendix <- data.frame(
  artifact = appendix_artifacts, scope = "structural_robustness",
  workflow_owner = ifelse(grepl("policy", appendix_artifacts), "optimal_policy",
    ifelse(grepl("campaign|randomization|reduced-form", appendix_artifacts),
           "reduced_form", "structural_postprocess")),
  default_contract = "frozen", source_path = appendix_artifacts,
  current_path = appendix_artifacts, stringsAsFactors = FALSE
)

registry <- rbind(
  main[c("artifact", "scope", "workflow_owner", "default_contract", "source_path", "current_path")],
  appendix
)
registry$distance_sensitive <- registry$default_contract == "generated" |
  registry$workflow_owner %in% c(
    "reduced_form", "balance", "structural_postprocess", "optimal_policy"
  )

candidate_roots <- c(
  file.path("build", "candidate-components", specification),
  file.path("build", "candidate-hpc", specification, "imported", "artifacts"),
  file.path("build", "structural-paper", specification, "fit105"),
  file.path("build", specification), file.path("build", "work", specification)
)
candidate_roots <- candidate_roots[dir.exists(candidate_roots)]
find_candidate <- function(artifact, source_path = NA_character_) {
  keys <- unique(na.omit(c(artifact, source_path)))
  keys <- sub("^[.]/", "", keys)
  hits <- unlist(lapply(candidate_roots, function(root) {
    paths <- list.files(root, recursive = TRUE, full.names = TRUE)
    relative <- substring(paths, nchar(root) + 2L)
    key <- sub("-1[.]pdf$", ".pdf", relative)
    exact <- paths[key %in% sub("-1[.]pdf$", ".pdf", keys)]
    if (length(exact)) return(exact)
    paths[sub("-1[.]pdf$", ".pdf", basename(paths)) %in%
            sub("-1[.]pdf$", ".pdf", basename(keys))]
  }), use.names = FALSE)
  unique(hits[file.exists(hits)])
}

dir.create(stage, recursive = TRUE, showWarnings = FALSE)
copy_one("ref-reports/ECM ReStud.tex", file.path(stage, "ECM ReStud.tex"))
appendix_destination <- file.path(stage, appendix_root)
dir.create(appendix_destination, recursive = TRUE, showWarnings = FALSE)
appendix_files <- list.files(appendix_root, recursive = TRUE, full.names = TRUE,
                             all.files = TRUE, no.. = TRUE)
appendix_files <- appendix_files[!dir.exists(appendix_files)]
for (source in appendix_files) {
  relative <- substring(source, nchar(appendix_root) + 2L)
  copy_one(source, file.path(appendix_destination, relative))
}
ri_section <- file.path(stage, appendix_root, "sections", "randomization-inference.tex")
if (file.exists(ri_section)) {
  ri_lines <- readLines(ri_section, warn = FALSE)
  ri_lines <- gsub(
    "six county by[[:space:]]*$", "six county by", ri_lines
  )
  ri_lines <- gsub("finalized-distance strata", "original assigned-distance strata", ri_lines,
                   fixed = TRUE)
  writeLines(ri_lines, ri_section)
}

registry$candidate_path <- NA_character_
registry$selected_path <- NA_character_
registry$selected_contract <- NA_character_
registry$status <- "missing"
for (index in seq_len(nrow(registry))) {
  candidates <- find_candidate(registry$artifact[[index]], registry$source_path[[index]])
  if (length(candidates) > 1L) {
    exact_suffix <- candidates[endsWith(candidates, registry$artifact[[index]])]
    if (length(exact_suffix) == 1L) candidates <- exact_suffix
  }
  if (length(candidates)) {
    source <- candidates[[1L]]
    contract_type <- "candidate"
    registry$candidate_path[[index]] <- source
  } else if (!registry$distance_sensitive[[index]] &&
             file.exists(registry$current_path[[index]])) {
    source <- registry$current_path[[index]]
    contract_type <- "verified_invariant_current"
  } else {
    next
  }
  destination <- file.path(stage, registry$artifact[[index]])
  if (copy_one(source, destination)) {
    registry$selected_path[[index]] <- source
    registry$selected_contract[[index]] <- contract_type
    registry$status[[index]] <- "staged"
  }
}

registry$current_sha256 <- vapply(registry$current_path, sha256, character(1))
registry$candidate_sha256 <- vapply(registry$candidate_path, sha256, character(1))
registry$changed <- !is.na(registry$candidate_sha256) &
  (is.na(registry$current_sha256) | registry$current_sha256 != registry$candidate_sha256)
dir.create(file.path(stage, "manifests"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(registry, file.path(stage, "artifact-comparison.csv"), row.names = FALSE, na = "")
missing <- registry$artifact[registry$status != "staged"]
writeLines(missing, file.path(stage, "missing-artifacts.txt"))
numeric_status <- system2("Rscript", c(
  "--no-save", "--no-restore", "scripts/checks/audit-candidate-numeric-prose.R",
  paste0("--stage=", stage)
))
if (!identical(numeric_status, 0L)) stop("Could not generate numerical prose TODOs.")
if (strict && length(missing)) {
  preliminary <- data.frame(
    distance_definition = specification,
    git_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]],
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    artifacts = nrow(registry), missing_artifacts = length(missing),
    manuscript_compiled = FALSE, appendix_compiled = FALSE,
    prose_status = "pending_numeric_todos"
  )
  utils::write.csv(preliminary, file.path(stage, "build-manifest.csv"), row.names = FALSE)
  jsonlite::write_json(as.list(preliminary[1, ]),
                       file.path(stage, "candidate-status.json"), auto_unbox = TRUE,
                       pretty = TRUE)
  stop("Candidate is missing ", length(missing), " required assigned artifacts. See ",
       file.path(stage, "missing-artifacts.txt"), call. = FALSE)
}

compile <- function(tex, workdir) {
  if (!nzchar(Sys.which("latexmk"))) return(FALSE)
  output <- system2("latexmk", c("-pdf", "-interaction=nonstopmode", basename(tex)),
                    stdout = TRUE, stderr = TRUE, env = character())
  status <- attr(output, "status")
  pdf <- sub("[.]tex$", ".pdf", basename(tex))
  ok <- is.null(status) || status == 0L
  if (!ok && file.exists(pdf)) {
    warning("LaTeX reported errors but produced ", file.path(workdir, pdf),
            "; inspect the log before promotion.")
    ok <- TRUE
  }
  ok
}
main_ok <- withr::with_dir(stage, compile("ECM ReStud.tex", stage))
appendix_stage <- file.path(stage, appendix_root)
appendix_ok <- withr::with_dir(appendix_stage,
  compile("structural-robustness-appendix.tex", appendix_stage))
appendix_pdf <- file.path(appendix_stage, "structural-robustness-appendix.pdf")
if (file.exists(appendix_pdf)) file.copy(appendix_pdf,
  file.path(stage, "structural-robustness-appendix.pdf"), overwrite = TRUE)

build_manifest <- data.frame(
  distance_definition = specification,
  git_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]],
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  artifacts = nrow(registry), missing_artifacts = length(missing),
  manuscript_compiled = main_ok, appendix_compiled = appendix_ok,
  prose_status = "pending_numeric_todos"
)
utils::write.csv(build_manifest, file.path(stage, "build-manifest.csv"), row.names = FALSE)
jsonlite::write_json(as.list(build_manifest[1, ]),
                     file.path(stage, "candidate-status.json"), auto_unbox = TRUE,
                     pretty = TRUE)
if (strict && (!main_ok || !appendix_ok)) stop("Candidate LaTeX compilation failed.")
cat(stage, "\n")
