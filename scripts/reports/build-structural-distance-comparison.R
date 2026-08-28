#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1L]])
}

fit_root <- value("--fit-root", "build/structural-fit/assigned")
output_root <- value("--output-root", "build/distance-comparison")
model <- "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
fit_pattern <- paste0("^dist_fit105_", model, "-[1-4][.]csv$")
chains <- sort(list.files(fit_root, pattern = fit_pattern, full.names = TRUE))
if (length(chains) != 4L) {
  stop("Expected four latest slim chain CSVs in ", fit_root, call. = FALSE)
}
chains <- normalizePath(chains)
chain_hashes <- vapply(chains, digest::digest, character(1L), file = TRUE,
                       algo = "sha256", serialize = FALSE)
chain_hashes <- unname(chain_hashes)
required_tables <- c(
  "struct-overall-te-table.tex", "private-signal-te-table.tex",
  "fob-beliefs-table.tex", "wtp-summ-table.tex"
)
manifest_path <- file.path(output_root, "structural-draw-manifest.csv")
expected <- unlist(lapply(c("realized", "assigned"), function(specification) {
  file.path(output_root, specification, "tables", required_tables)
}))

if (file.exists(manifest_path) && all(file.exists(expected))) {
  prior <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (identical(prior$chain, chains) &&
      identical(prior$sha256, chain_hashes)) {
    message("Structural comparison outputs already match the current chain hashes.")
    quit(save = "no", status = 0L)
  }
}

run <- function(script, script_args, log_path) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript", c("--no-save", "--no-restore", script, script_args),
    stdout = log_path, stderr = log_path
  )
  if (!identical(status, 0L)) {
    stop("Command failed; see ", log_path, call. = FALSE)
  }
}

for (specification in c("realized", "assigned")) {
  root <- file.path(output_root, specification)
  compact <- file.path(root, "structural", "compact-gq")
  data <- file.path(root, "structural-data")
  tables <- file.path(root, "tables")
  figures <- file.path(root, "figures")
  dir.create(compact, recursive = TRUE, showWarnings = FALSE)
  dir.create(data, recursive = TRUE, showWarnings = FALSE)

  run(
    "scripts/structural/generate-compact-gq.R",
    c(
      "--workspace=build/structural-workspace/main-core-input.RData",
      paste0("--fit-csvs=", paste(chains, collapse = ",")),
      paste0("--output-path=", compact),
      paste0("--distance-definition=", specification)
    ),
    file.path(root, "compact-gq.log")
  )
  run(
    "scripts/structural/postprocess-main-core-compact.R",
    c(
      paste0("--compact-gq-path=", compact),
      paste0("--fit-path=", fit_root),
      paste0("--output-path=", data),
      "--fit-version=105",
      paste0("--model=", model)
    ),
    file.path(root, "structural-postprocess.log")
  )
  run(
    "scripts/structural/render-paper.R",
    c(
      "--fit-version=105",
      paste0("--model=", model),
      paste0("--input-path=", data),
      paste0("--output-path=", data),
      paste0("--table-output=", tables),
      paste0("--figure-output=", figures),
      "--tables-only"
    ),
    file.path(root, "structural-render.log")
  )
  missing <- file.path(tables, required_tables)
  missing <- missing[!file.exists(missing)]
  if (length(missing)) stop("Missing structural comparison tables: ",
                            paste(missing, collapse = ", "), call. = FALSE)
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(
    chain = chains,
    sha256 = chain_hashes,
    posterior_role = "shared by realized and assigned postprocessing",
    stringsAsFactors = FALSE
  ),
  manifest_path,
  row.names = FALSE
)
message("Wrote structural distance comparison outputs to ", output_root)
