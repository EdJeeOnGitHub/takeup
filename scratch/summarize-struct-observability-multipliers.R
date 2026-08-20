#!/usr/bin/env Rscript

# Compare the baseline social multiplier with completed structural fits that
# replace the first-order definite-classification observability input.  The
# missing-data-marginalized row is analytically identical to baseline in the
# current observability model: missing rows have no likelihood contribution,
# and the predicted probability depends only on arm and cluster distance, so
# averaging predicted probabilities over observed plus missing respondents
# returns the same cluster schedule.

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Option supplied more than once: ", name,
                             call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

input_path <- option_value("--input-path", "temp-data/struct-postprocess")
csv_output <- option_value(
  "--csv-output",
  "temp-data/struct-observability-multiplier-summary.csv"
)
tex_output <- option_value(
  "--tex-output",
  paste0("appendix/structural-robustness/tables/",
         "main-core-observability-multipliers.tex")
)

specifications <- data.frame(
  fit = c(105L, 106L, 106L),
  model = c(
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP",
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS",
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB"
  ),
  specification = c(
    "First-order definite classification (baseline)",
    "Correct classification",
    "Perceived observability (second order)"
  ),
  stringsAsFactors = FALSE
)

summary_path <- function(fit, model) file.path(
  input_path,
  sprintf("rvar_processed_dist_fit%d_sm_summ_%s_1-4.rds", fit, model)
)
paths <- mapply(summary_path, specifications$fit, specifications$model,
                USE.NAMES = FALSE)
if (any(!file.exists(paths))) {
  stop("Missing multiplier summaries: ",
       paste(basename(paths[!file.exists(paths)]), collapse = ", "),
       call. = FALSE)
}

treatments <- c("Bracelet", "Calendar", "Ink", "Control")
distances <- c(500, 1500, 2500)
rows <- vector("list", nrow(specifications))
for (i in seq_len(nrow(specifications))) {
  x <- readRDS(paths[[i]])
  required <- c("model", "fit_version", "treatment", "roc_distance",
                "variable", "value", "conf.low", "conf.high", ".width")
  if (!all(required %in% names(x))) {
    stop("Malformed summary: ", paths[[i]], call. = FALSE)
  }
  if (!identical(unique(as.character(x$model)), specifications$model[[i]]) ||
      !identical(unique(as.integer(x$fit_version)), specifications$fit[[i]])) {
    stop("Model or fit mismatch in: ", paths[[i]], call. = FALSE)
  }
  x <- x[
    x$variable == "sm_rescaled" &
      x$roc_distance %in% distances &
      as.character(x$treatment) %in% treatments,
    , drop = FALSE
  ]
  if (nrow(x) != length(treatments) * length(distances) ||
      any(x$.width != 0.95)) {
    stop("Incomplete 95% multiplier grid in: ", paths[[i]], call. = FALSE)
  }
  rows[[i]] <- data.frame(
    specification = specifications$specification[[i]],
    treatment = as.character(x$treatment),
    distance = x$roc_distance,
    median = -x$value,
    lower = -x$conf.high,
    upper = -x$conf.low,
    source_fit = specifications$fit[[i]],
    source_model = specifications$model[[i]],
    stringsAsFactors = FALSE
  )
}
results <- do.call(rbind, rows)

# Under the current specification this is not a fourth fitted posterior. It is
# the baseline target expressed with explicit marginalization of belief rows
# whose outcomes are missing. Retain that fact in the machine-readable output.
missing_rows <- results[
  results$specification == specifications$specification[[1L]],
  , drop = FALSE
]
missing_rows$specification <- "First-order observability, missing data marginalized"
missing_rows$source_model <- "Analytically identical to baseline target"
results <- rbind(results, missing_rows)

specification_order <- c(
  specifications$specification,
  "First-order observability, missing data marginalized"
)
results$specification <- factor(results$specification,
                                levels = specification_order)
results$treatment <- factor(results$treatment, levels = treatments)
results <- results[order(results$specification, results$treatment,
                         results$distance), ]
results$specification <- as.character(results$specification)
results$treatment <- as.character(results$treatment)

dir.create(dirname(csv_output), recursive = TRUE, showWarnings = FALSE)
write.csv(results, csv_output, row.names = FALSE)

fmt <- function(x) sprintf("%.2f", x)
ci <- function(lo, hi) paste0("(", fmt(lo), ", ", fmt(hi), ")")
panel_titles <- c(
  "Panel A: First-order definite classification (baseline)",
  "Panel B: Correct classification",
  "Panel C: Perceived observability (second order)",
  "Panel D: First-order observability, missing data marginalized"
)
lines <- c(
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Treatment & 500m & 1500m & 2500m \\\\",
  "\\midrule"
)
for (i in seq_along(specification_order)) {
  block <- results[results$specification == specification_order[[i]], ]
  if (i > 1L) lines <- c(lines, "\\addlinespace")
  lines <- c(lines, paste0(
    "\\multicolumn{4}{l}{\\textit{", panel_titles[[i]], "}} \\\\"
  ))
  for (treatment_name in treatments) {
    cells <- block[block$treatment == treatment_name, ]
    cells <- cells[match(distances, cells$distance), ]
    lines <- c(
      lines,
      paste0(treatment_name, " & ",
             paste(fmt(cells$median), collapse = " & "), " \\\\"),
      paste0(" & ", paste(ci(cells$lower, cells$upper), collapse = " & "),
             " \\\\ ")
    )
  }
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
dir.create(dirname(tex_output), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, tex_output)

message("Wrote: ", csv_output)
message("Wrote: ", tex_output)
