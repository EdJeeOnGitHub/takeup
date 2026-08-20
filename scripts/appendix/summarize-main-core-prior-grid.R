#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")

manifest_path <- main_core_option_value(args, "--manifest")
output_path <- main_core_option_value(
  args, "--output-path", "temp-data/main-core-prior-grid"
)
table_path <- main_core_option_value(
  args, "--table-path", file.path(output_path, "main-core-prior-grid.tex")
)
if (is.null(manifest_path) || !file.exists(manifest_path)) {
  stop("--manifest must name an existing manifest.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(manifest) != 13L || anyDuplicated(manifest$label)) {
  stop("Expected the unique 13-specification prior manifest.", call. = FALSE)
}
treatments <- c("Control", "Ink", "Calendar", "Bracelet")
spans <- c("Combined", "Short", "Long")
all_draws <- list()
all_summaries <- list()
all_diagnostics <- list()
all_ates <- list()

summarize_value <- function(value) c(
  mean = mean(value), sd = sd(value),
  q025 = unname(quantile(value, .025)),
  q50 = unname(quantile(value, .5)),
  q975 = unname(quantile(value, .975))
)

for (spec_index in seq_len(nrow(manifest))) {
  specification <- manifest[spec_index, ]
  run_candidates <- data.frame(
    fit_root = c("fits-rerun-td15", "fits-rerun", "fits"),
    gq_root = c("gq-rerun-td15", "gq-rerun", "gq"),
    run_source = c(
      "800/800 treedepth-15 rerun", "800/800 rerun", "400/400 initial"
    ),
    max_treedepth_limit = c(15L, 12L, 12L),
    stringsAsFactors = FALSE
  )
  complete_candidate <- vapply(seq_len(nrow(run_candidates)), function(i) {
    length(list.files(
      file.path(output_path, run_candidates$fit_root[i], specification$label),
      pattern = "\\.csv$", full.names = TRUE
    )) == 4L && length(list.files(
      file.path(output_path, run_candidates$gq_root[i], specification$label),
      pattern = "\\.csv$", full.names = TRUE
    )) == 4L
  }, logical(1))
  if (!any(complete_candidate)) {
    stop("No complete fit/GQ pair for ", specification$label, call. = FALSE)
  }
  run <- run_candidates[which(complete_candidate)[1L], ]
  fit_dir <- file.path(output_path, run$fit_root, specification$label)
  gq_dir <- file.path(output_path, run$gq_root, specification$label)
  fit_files <- sort(list.files(
    fit_dir, pattern = "\\.csv$", full.names = TRUE
  ))
  gq_files <- sort(list.files(
    gq_dir, pattern = "\\.csv$", full.names = TRUE
  ))
  if (length(fit_files) != 4L || length(gq_files) != 4L) {
    stop("Expected four fit and four GQ CSVs for ", specification$label,
         call. = FALSE)
  }
  fit <- read_cmdstan_csv(fit_files)
  gq_fit <- read_cmdstan_csv(gq_files)
  fit_summary <- summarise_draws(fit$post_warmup_draws)
  sampler <- as_draws_df(fit$post_warmup_sampler_diagnostics)
  gq <- as_draws_df(gq_fit$generated_quantities)
  finite_names <- as.vector(outer(
    seq_len(3L), seq_len(4L),
    function(span, treatment) sprintf(
      "core_compact_gaussian_sm_finite_segments[%d,%d]", span, treatment
    )
  ))
  if (!all(finite_names %in% names(gq))) {
    stop("Finite-segment GQ is incomplete for ", specification$label,
         call. = FALSE)
  }
  finite_diagnostics <- summarise_draws(subset_draws(
    gq_fit$generated_quantities,
    variable = finite_names
  ))
  all_diagnostics[[spec_index]] <- data.frame(
    spec_id = specification$spec_id,
    label = specification$label,
    panel = specification$panel,
    setting = specification$setting,
    run_source = run$run_source,
    max_rhat_parameters = max(fit_summary$rhat, na.rm = TRUE),
    min_ess_bulk_parameters = min(fit_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail_parameters = min(fit_summary$ess_tail, na.rm = TRUE),
    max_rhat_estimands = max(finite_diagnostics$rhat, na.rm = TRUE),
    min_ess_bulk_estimands = min(finite_diagnostics$ess_bulk, na.rm = TRUE),
    min_ess_tail_estimands = min(finite_diagnostics$ess_tail, na.rm = TRUE),
    divergences = sum(sampler$divergent__, na.rm = TRUE),
    max_treedepth = sum(
      sampler$treedepth__ >= run$max_treedepth_limit, na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )

  draw_id <- seq_len(nrow(gq))
  draw_rows <- list()
  for (span_index in seq_along(spans)) {
    values <- sapply(seq_along(treatments), function(treatment_index) {
      gq[[sprintf(
        "core_compact_gaussian_sm_finite_segments[%d,%d]",
        span_index, treatment_index
      )]]
    })
    colnames(values) <- treatments
    no_signal <- rowMeans(values[, c("Control", "Calendar"), drop = FALSE])
    signal <- rowMeans(values[, c("Ink", "Bracelet"), drop = FALSE])
    arm_rows <- do.call(rbind, lapply(seq_along(treatments), function(index) {
      data.frame(
        spec_id = specification$spec_id,
        label = specification$label,
        panel = specification$panel,
        setting = specification$setting,
        draw = draw_id,
        span = spans[span_index],
        estimand = "Arm multiplier",
        treatment = treatments[index],
        value = values[, index],
        stringsAsFactors = FALSE
      )
    }))
    comparison_rows <- rbind(
      data.frame(
        spec_id = specification$spec_id, label = specification$label,
        panel = specification$panel, setting = specification$setting,
        draw = draw_id, span = spans[span_index],
        estimand = "No Signal multiplier", treatment = "No Signal",
        value = no_signal
      ),
      data.frame(
        spec_id = specification$spec_id, label = specification$label,
        panel = specification$panel, setting = specification$setting,
        draw = draw_id, span = spans[span_index],
        estimand = "Signal multiplier", treatment = "Signal", value = signal
      ),
      data.frame(
        spec_id = specification$spec_id, label = specification$label,
        panel = specification$panel, setting = specification$setting,
        draw = draw_id, span = spans[span_index],
        estimand = "No Signal - Signal", treatment = "Comparison",
        value = no_signal - signal
      ),
      data.frame(
        spec_id = specification$spec_id, label = specification$label,
        panel = specification$panel, setting = specification$setting,
        draw = draw_id, span = spans[span_index],
        estimand = "Calendar - Bracelet", treatment = "Comparison",
        value = values[, "Calendar"] - values[, "Bracelet"]
      )
    )
    draw_rows[[span_index]] <- rbind(arm_rows, comparison_rows)
  }
  specification_draws <- do.call(rbind, draw_rows)
  all_draws[[spec_index]] <- specification_draws

  summary_groups <- split(
    specification_draws,
    interaction(
      specification_draws$span, specification_draws$estimand,
      specification_draws$treatment, drop = TRUE
    )
  )
  all_summaries[[spec_index]] <- do.call(rbind, lapply(summary_groups, function(x) {
    summary <- summarize_value(x$value)
    data.frame(
      spec_id = x$spec_id[1], label = x$label[1], panel = x$panel[1],
      setting = x$setting[1], span = x$span[1],
      estimand = x$estimand[1], treatment = x$treatment[1],
      t(summary),
      directional_probability = if (x$estimand[1] %in% c(
        "No Signal - Signal", "Calendar - Bracelet"
      )) mean(x$value > 0) else NA_real_,
      n_draws = nrow(x), check.names = FALSE
    )
  }))

  takeup <- sapply(seq_along(treatments), function(treatment_index) {
    gq[[sprintf("core_compact_takeup_level[1,%d]", treatment_index)]]
  })
  colnames(takeup) <- treatments
  for (treatment in c("Ink", "Calendar", "Bracelet")) {
    value <- takeup[, treatment] - takeup[, "Control"]
    all_ates[[length(all_ates) + 1L]] <- data.frame(
      spec_id = specification$spec_id, label = specification$label,
      treatment = treatment, t(summarize_value(value)), check.names = FALSE
    )
  }
}

draws <- do.call(rbind, all_draws)
summary <- do.call(rbind, all_summaries)
diagnostics <- do.call(rbind, all_diagnostics)
ates <- do.call(rbind, all_ates)
diagnostics$needs_rerun <- with(
  diagnostics,
  max_rhat_parameters > 1.01 | max_rhat_estimands > 1.01 |
    min_ess_bulk_parameters < 400 | min_ess_tail_parameters < 400 |
    min_ess_bulk_estimands < 400 | min_ess_tail_estimands < 400 |
    divergences > 0 | max_treedepth > 0
)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
write.csv(draws, file.path(output_path, "prior-grid-estimand-draws.csv"),
          row.names = FALSE)
write.csv(summary, file.path(output_path, "prior-grid-estimand-summary.csv"),
          row.names = FALSE)
write.csv(ates, file.path(output_path, "prior-grid-ate-summary.csv"),
          row.names = FALSE)
write.csv(diagnostics, file.path(output_path, "prior-grid-diagnostics.csv"),
          row.names = FALSE)
write.csv(diagnostics[diagnostics$needs_rerun, ],
          file.path(output_path, "prior-grid-needs-rerun.csv"), row.names = FALSE)

paper <- summary[summary$span == "Combined", ]
cell <- function(label, treatment) {
  x <- paper[
    paper$label == label & paper$estimand == "Arm multiplier" &
      paper$treatment == treatment,
  ]
  if (nrow(x) != 1L) stop("Missing table cell: ", label, "/", treatment)
  sprintf("\\makecell[c]{%.2f\\\\(%.2f, %.2f)}", x$q50, x$q025, x$q975)
}
probability_cell <- function(label, estimand) {
  x <- paper[paper$label == label & paper$estimand == estimand, ]
  if (nrow(x) != 1L) stop("Missing probability cell: ", label, "/", estimand)
  if (x$directional_probability >= 0.9995) ">0.999" else
    sprintf("%.3f", x$directional_probability)
}

lines <- c(
  "\\centering",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
  "\\begin{tabular}[t]{lcccccc}", "\\toprule",
  paste0(
    " & \\multicolumn{4}{c}{Finite multiplier, 0--2.5 km} & ",
    "\\multicolumn{2}{c}{Directional posterior probability} \\\\"
  ),
  "\\cmidrule(l{3pt}r{3pt}){2-5}\\cmidrule(l{3pt}r{3pt}){6-7}",
  paste0(
    "Prior specification & Bracelet & Calendar & Ink & Control & ",
    "$M_{NS}>M_S$ & $M_C>M_B$ \\\\"
  ),
  "\\midrule"
)
panel_order <- unique(manifest$panel)
for (panel_index in seq_along(panel_order)) {
  panel <- panel_order[panel_index]
  block <- manifest[manifest$panel == panel, ]
  if (panel != "Baseline") {
    lines <- c(
      lines,
      paste0("\\multicolumn{7}{l}{\\textit{", panel, "}}\\\\"),
      "\\addlinespace"
    )
  }
  for (row_index in seq_len(nrow(block))) {
    specification <- block[row_index, ]
    lines <- c(lines, paste0(
      specification$setting, " & ",
      paste(vapply(
        c("Bracelet", "Calendar", "Ink", "Control"),
        function(treatment) cell(specification$label, treatment), ""
      ), collapse = " & "), " & ",
      probability_cell(specification$label, "No Signal - Signal"), " & ",
      probability_cell(specification$label, "Calendar - Bracelet"),
      "\\\\"
    ))
    if (row_index < nrow(block)) lines <- c(lines, "\\addlinespace")
  }
  if (panel_index < length(panel_order)) lines <- c(lines, "\\addlinespace")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}%", "}")
writeLines(lines, table_path)

message("Wrote prior-grid table: ", table_path)
message("Specifications needing longer reruns: ", sum(diagnostics$needs_rerun))
