#!/usr/bin/env Rscript

# Build the structural social-multiplier robustness table from the compact
# summaries produced by scripts/structural/postprocess-roc.R --sm. The corresponding draw
# objects and all four raw chains are validated before any paper-facing output
# is written.
#
# Pass --individual-fp-production-gq-path=DIR to replace the legacy short-run
# individual-FP row with all compact production GQ chain files found in DIR.
# This supports a provisional two-chain table and needs no code change when DIR
# is later regenerated with all four chains.

args <- commandArgs(trailingOnly = TRUE)

option_value <- function(name, default) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) {
    stop("Option supplied more than once: ", name, call. = FALSE)
  }
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

has_flag <- function(name) name %in% args

input_path <- option_value("--input-path", "temp-data/struct-postprocess")
raw_path <- option_value("--raw-path", "data/stan_analysis_data")
compact_gq_path <- option_value(
  "--compact-gq-path",
  file.path(input_path, "compact-sm-gq")
)
individual_fp_production_gq_path <- option_value(
  "--individual-fp-production-gq-path",
  ""
)
convention_path <- option_value(
  "--convention-path",
  "temp-data/social-multiplier-decomposition-values.csv"
)
csv_output <- option_value(
  "--csv-output",
  "temp-data/struct-social-multiplier-robustness-summary.csv"
)
tex_output <- option_value(
  "--tex-output",
  "presentations/tables/fit105/struct-social-multiplier-robustness-table.tex"
)
include_secondary <- has_flag("--include-secondary")
allow_incomplete <- has_flag("--allow-incomplete")
skip_draw_validation <- has_flag("--skip-draw-validation")

unknown <- args[
  !startsWith(args, "--input-path=") &
    !startsWith(args, "--raw-path=") &
    !startsWith(args, "--compact-gq-path=") &
    !startsWith(args, "--individual-fp-production-gq-path=") &
    !startsWith(args, "--convention-path=") &
    !startsWith(args, "--csv-output=") &
    !startsWith(args, "--tex-output=") &
    !args %in% c(
      "--include-secondary",
      "--allow-incomplete",
      "--skip-draw-validation"
    )
]
if (length(unknown) > 0L) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

models <- data.frame(
  fit = c(105L, 105L, 105L, 105L, 106L, 106L),
  model = c(
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP",
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS",
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP",
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS",
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB"
  ),
  specification = c(
    "Baseline Specification",
    "Private Distance Costs and Community Social Image Returns",
    "Individual Distance Observed by Peers",
    "Excluding Spatially Dispersed Clusters",
    "Corrected observability",
    "SOB observability"
  ),
  secondary = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  compact_gq = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
models <- models[!models$secondary | include_secondary, , drop = FALSE]

required_quantities <- c(
  "cluster_social_multiplier",
  "cluster_mu_rep",
  "cluster_mu_rep_deriv",
  "cluster_delta",
  "cluster_delta_deriv",
  "dist_beta_v"
)
required_treatments <- c("Bracelet", "Calendar", "Ink", "Control")
target_distances <- c(500, 1500, 2500)
required_variables <- c(
  "sm",
  "dist_beta_v",
  "sm_delta_part",
  "sm_mu_part",
  "sm_rescaled",
  "sm_delta_part_rescaled",
  "sm_mu_part_rescaled"
)

object_path <- function(fit, model, kind) {
  file.path(
    input_path,
    sprintf(
      "rvar_processed_dist_fit%d_sm_%s_%s_1-4.rds",
      fit,
      kind,
      model
    )
  )
}

chain_path <- function(chain, fit, model) {
  file.path(
    raw_path,
    sprintf("dist_fit%d_%s-%d.csv", fit, model, chain)
  )
}

compact_chain_path <- function(chain, model) {
  file.path(
    compact_gq_path,
    sprintf("compact_sm_fit105_%s-%d.csv", model, chain)
  )
}

read_stan_header <- function(path) {
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      stop("No Stan CSV header found in: ", path, call. = FALSE)
    }
    if (!startsWith(line, "#")) return(strsplit(line, ",", fixed = TRUE)[[1L]])
  }
}

base_quantity <- function(x) sub("\\..*$", "", x)

preflight_model <- function(fit, model, compact_gq) {
  chains <- vapply(1:4, chain_path, character(1), fit = fit, model = model)
  missing_chains <- chains[!file.exists(chains)]
  if (length(missing_chains) > 0L) {
    return(paste("missing chain(s):", paste(basename(missing_chains), collapse = ", ")))
  }

  if (compact_gq) {
    compact_chains <- vapply(
      1:4,
      compact_chain_path,
      character(1),
      model = model
    )
    missing_compact <- compact_chains[!file.exists(compact_chains)]
    if (length(missing_compact) > 0L) {
      return(paste(
        "missing compact GQ chain(s):",
        paste(basename(missing_compact), collapse = ", ")
      ))
    }
    compact_variables <- paste0("compact_", required_variables)
    compact_variables[compact_variables == "compact_dist_beta_v"] <-
      "compact_dist_beta_v"
    for (chain in seq_along(compact_chains)) {
      header_quantities <- unique(base_quantity(
        read_stan_header(compact_chains[[chain]])
      ))
      missing_quantities <- setdiff(compact_variables, header_quantities)
      if (length(missing_quantities) > 0L) {
        return(sprintf(
          "compact GQ chain %d lacks quantities: %s",
          chain,
          paste(missing_quantities, collapse = ", ")
        ))
      }
    }
  } else {
    for (chain in seq_along(chains)) {
      header_quantities <- unique(base_quantity(read_stan_header(chains[[chain]])))
      missing_quantities <- setdiff(required_quantities, header_quantities)
      if (length(missing_quantities) > 0L) {
        return(sprintf(
          "chain %d lacks generated quantities: %s",
          chain,
          paste(missing_quantities, collapse = ", ")
        ))
      }
    }
  }

  missing_objects <- c(
    object_path(fit, model, "summ"),
    object_path(fit, model, "draws")
  )
  missing_objects <- missing_objects[!file.exists(missing_objects)]
  if (length(missing_objects) > 0L) {
    return(paste(
      "missing scripts/structural/postprocess-roc.R --sm output(s):",
      paste(basename(missing_objects), collapse = ", ")
    ))
  }
  NA_character_
}

preflight <- vapply(
  seq_len(nrow(models)),
  function(i) preflight_model(
    models$fit[[i]],
    models$model[[i]],
    models$compact_gq[[i]]
  ),
  character(1)
)

if (any(!is.na(preflight))) {
  failures <- sprintf(
    "- %s (fit %d): %s",
    models$specification[!is.na(preflight)],
    models$fit[!is.na(preflight)],
    preflight[!is.na(preflight)]
  )
  message(
    "Unavailable structural social-multiplier specifications:\n",
    paste(failures, collapse = "\n")
  )
  if (!allow_incomplete) {
    stop(
      "No paper-facing outputs were written. Supply the missing generated ",
      "quantities/postprocessed objects, or use --allow-incomplete only for ",
      "a diagnostic partial table.",
      call. = FALSE
    )
  }
  models <- models[is.na(preflight), , drop = FALSE]
  if (identical(
    csv_output,
    "temp-data/struct-social-multiplier-robustness-summary.csv"
  )) {
    csv_output <- "temp-data/struct-social-multiplier-robustness-summary-incomplete.csv"
  }
  if (identical(
    tex_output,
    "presentations/tables/fit105/struct-social-multiplier-robustness-table.tex"
  )) {
    tex_output <- paste0(
      "presentations/tables/fit105/",
      "struct-social-multiplier-robustness-table-incomplete.tex"
    )
  }
}

if (nrow(models) == 0L) {
  stop("No complete model outputs are available.", call. = FALSE)
}

validate_grid <- function(x, label, expected_fit, expected_model) {
  needed_columns <- c(
    "model", "fit_version", "fit_type", "treatment", "roc_distance", "variable"
  )
  missing_columns <- setdiff(needed_columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(label, " lacks columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (!identical(unique(as.integer(x$fit_version)), as.integer(expected_fit))) {
    stop(label, " has the wrong fit_version.", call. = FALSE)
  }
  if (!identical(unique(as.character(x$model)), expected_model)) {
    stop(label, " has the wrong model name.", call. = FALSE)
  }
  if (!identical(unique(as.character(x$fit_type)), "fit")) {
    stop(label, " is not a posterior-fit object.", call. = FALSE)
  }

  treatments <- unique(as.character(x$treatment))
  if (!setequal(treatments, required_treatments)) {
    stop(
      label, " treatment labels differ from the required four: ",
      paste(treatments, collapse = ", "),
      call. = FALSE
    )
  }
  variables <- unique(as.character(x$variable))
  if (!setequal(variables, required_variables)) {
    stop(
      label, " variables differ from the expected social-multiplier set: ",
      paste(variables, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(x$roc_distance > 5000 | x$roc_distance < 0)) {
    stop(label, " contains distances outside [0, 5000].", call. = FALSE)
  }

  keys <- paste(x$treatment, x$roc_distance, x$variable, sep = "\r")
  if (anyDuplicated(keys)) {
    stop(label, " has duplicate treatment-distance-variable rows.", call. = FALSE)
  }
  expected_rows <- length(required_treatments) *
    length(unique(x$roc_distance)) *
    length(required_variables)
  if (nrow(x) != expected_rows) {
    stop(
      label, " has ", nrow(x), " rows; expected ", expected_rows,
      " for a complete treatment-distance-variable grid.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

summaries <- vector("list", nrow(models))
for (i in seq_len(nrow(models))) {
  model_label <- sprintf("%s (fit %d)", models$specification[[i]], models$fit[[i]])
  summary_path <- object_path(models$fit[[i]], models$model[[i]], "summ")
  summary_object <- readRDS(summary_path)
  validate_grid(
    summary_object,
    paste(model_label, "summary"),
    models$fit[[i]],
    models$model[[i]]
  )
  summary_columns <- c("value", "conf.low", "conf.high", ".width")
  if (!all(summary_columns %in% names(summary_object))) {
    stop(
      model_label, " summary lacks interval columns: ",
      paste(setdiff(summary_columns, names(summary_object)), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(summary_object$.width != 0.95)) {
    stop(model_label, " summary does not use 95% credible intervals.", call. = FALSE)
  }
  summaries[[i]] <- summary_object

  if (!skip_draw_validation) {
    draw_path <- object_path(models$fit[[i]], models$model[[i]], "draws")
    message("Reading and validating draw object: ", draw_path)
    draw_object <- readRDS(draw_path)
    validate_grid(
      draw_object,
      paste(model_label, "draws"),
      models$fit[[i]],
      models$model[[i]]
    )
    rm(draw_object)
    invisible(gc())
  }
}

# Optionally replace the legacy short-run individual-fixed-point multiplier
# rows with a compact GQ run from the production fit. This deliberately reads
# however many compact chain files are present: two can be used provisionally,
# and the identical command automatically uses all four after the final GQ run.
if (nzchar(individual_fp_production_gq_path)) {
  individual_fp_model <-
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP"
  individual_fp_index <- which(models$model == individual_fp_model)
  if (length(individual_fp_index) != 1L) {
    stop("The individual fixed-point specification is unavailable.", call. = FALSE)
  }

  production_files <- list.files(
    individual_fp_production_gq_path,
    pattern = paste0(
      "^compact_sm_fit105_", individual_fp_model, "-[0-9]+[.]csv$"
    ),
    full.names = TRUE
  )
  if (length(production_files) < 2L) {
    stop(
      "Expected at least two production individual-FP compact GQ files in ",
      individual_fp_production_gq_path,
      call. = FALSE
    )
  }
  production_draws <- do.call(rbind, lapply(production_files, function(path) {
    read.csv(path, comment.char = "#", check.names = FALSE)
  }))
  stan_treatments <- c("Control", "Ink", "Calendar", "Bracelet")
  replacement <- summaries[[individual_fp_index]]
  valid_counts <- integer(0)
  for (distance_index in seq_along(target_distances)) {
    for (treatment_index in seq_along(stan_treatments)) {
      column <- sprintf(
        "compact_sm_rescaled.%d.%d",
        distance_index,
        treatment_index
      )
      if (!column %in% names(production_draws)) {
        stop("Production compact GQ lacks column: ", column, call. = FALSE)
      }
      values <- production_draws[[column]]
      values <- values[is.finite(values)]
      if (length(values) == 0L) {
        stop("Production compact GQ has no finite draws for: ", column, call. = FALSE)
      }
      valid_counts <- c(valid_counts, length(values))
      row <- replacement$variable == "sm_rescaled" &
        replacement$roc_distance == target_distances[[distance_index]] &
        as.character(replacement$treatment) == stan_treatments[[treatment_index]]
      if (sum(row) != 1L) {
        stop("Could not uniquely locate the individual-FP summary row.", call. = FALSE)
      }
      replacement$value[row] <- median(values)
      replacement$conf.low[row] <- unname(quantile(values, 0.025))
      replacement$conf.high[row] <- unname(quantile(values, 0.975))
    }
  }
  summaries[[individual_fp_index]] <- replacement
  message(
    "Using production individual-FP compact GQ: ",
    length(production_files), " chain file(s), ", min(valid_counts),
    " finite draws per reported cell."
  )
}

main_index <- which(models$model == "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP")
if (length(main_index) != 1L) {
  stop("The main structural model is required to establish the sign convention.", call. = FALSE)
}
if (!file.exists(convention_path)) {
  stop("Existing paper-convention file not found: ", convention_path, call. = FALSE)
}

paper_values <- read.csv(convention_path, check.names = FALSE)
paper_column <- "Social Multiplier"
paper_required <- c("treatment", "distance", paper_column)
if (!all(paper_required %in% names(paper_values))) {
  stop(
    "Convention file lacks columns: ",
    paste(setdiff(paper_required, names(paper_values)), collapse = ", "),
    call. = FALSE
  )
}

main_summary <- summaries[[main_index]]
main_sm <- main_summary[
  main_summary$variable == "sm_rescaled" &
    main_summary$roc_distance <= 2500,
  c("treatment", "roc_distance", "value")
]
names(main_sm)[names(main_sm) == "roc_distance"] <- "distance"
comparison <- merge(
  main_sm,
  paper_values[, paper_required],
  by = c("treatment", "distance")
)
if (nrow(comparison) == 0L) {
  stop("No overlapping rows were available to establish the sign convention.", call. = FALSE)
}

same_error <- max(abs(comparison$value - comparison[[paper_column]]))
flipped_error <- max(abs(-comparison$value - comparison[[paper_column]]))
sign_multiplier <- if (same_error <= flipped_error) 1 else -1
convention_error <- min(same_error, flipped_error)
if (!is.finite(convention_error) || convention_error > 1e-8) {
  stop(
    "Main-model sm_rescaled does not reproduce the existing paper convention; ",
    "maximum absolute discrepancy is ", signif(convention_error, 4), ".",
    call. = FALSE
  )
}
message(
  "Paper sign convention confirmed: reported multiplier = ",
  if (sign_multiplier == -1) "-sm_rescaled" else "sm_rescaled",
  "."
)

nearest_distance <- function(target, available) {
  available[[which.min(abs(available - target))]]
}

selected_rows <- vector("list", nrow(models))
for (i in seq_len(nrow(models))) {
  x <- summaries[[i]]
  x <- x[
    x$variable == "sm_rescaled" &
      x$roc_distance <= 2500 &
      as.character(x$treatment) %in% required_treatments,
    ,
    drop = FALSE
  ]
  available <- sort(unique(x$roc_distance))
  selected <- vapply(target_distances, nearest_distance, numeric(1), available = available)
  if (any(abs(selected - target_distances) > 100)) {
    stop(
      models$specification[[i]],
      " has no distance within 100m of one or more requested distances.",
      call. = FALSE
    )
  }

  x <- x[x$roc_distance %in% selected, , drop = FALSE]
  x$target_distance <- target_distances[match(x$roc_distance, selected)]
  x$Specification <- models$specification[[i]]
  x$Treatment <- as.character(x$treatment)
  x$median <- sign_multiplier * x$value
  if (sign_multiplier == 1) {
    x$lower <- x$conf.low
    x$upper <- x$conf.high
  } else {
    x$lower <- -x$conf.high
    x$upper <- -x$conf.low
  }
  selected_rows[[i]] <- x
}
selected <- do.call(rbind, selected_rows)
expected_selected_rows <- nrow(models) *
  length(required_treatments) *
  length(target_distances)
if (nrow(selected) != expected_selected_rows) {
  stop(
    "Selected multiplier grid has ", nrow(selected), " rows; expected ",
    expected_selected_rows, ".", call. = FALSE
  )
}

format_number <- function(x) {
  out <- formatC(x, format = "f", digits = 3)
  sub("-0\\.000", "0.000", out)
}
selected$cell <- sprintf(
  "%s [%s, %s]",
  format_number(selected$median),
  format_number(selected$lower),
  format_number(selected$upper)
)

row_order <- expand.grid(
  Treatment = required_treatments,
  Specification = models$specification,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
wide <- reshape(
  selected[, c("Specification", "Treatment", "target_distance", "cell")],
  idvar = c("Specification", "Treatment"),
  timevar = "target_distance",
  direction = "wide"
)
names(wide) <- sub("^cell\\.", "M(", names(wide))
names(wide) <- sub("^(M\\([0-9]+)$", "\\1m)", names(wide))
wide <- merge(row_order, wide, by = c("Specification", "Treatment"), sort = FALSE)
wide <- wide[
  match(
    paste(row_order$Specification, row_order$Treatment, sep = "\r"),
    paste(wide$Specification, wide$Treatment, sep = "\r")
  ),
  c("Specification", "Treatment", "M(500m)", "M(1500m)", "M(2500m)")
]

dir.create(dirname(csv_output), recursive = TRUE, showWarnings = FALSE)
write.csv(wide, csv_output, row.names = FALSE, na = "")

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)
  x
}
panel_letters <- LETTERS[seq_len(nrow(models))]
latex_rows <- character(0)
for (model_index in seq_len(nrow(models))) {
  specification <- models$specification[[model_index]]
  if (model_index > 1L) {
    latex_rows <- c(latex_rows, "\\addlinespace[0.3em]")
  }
  panel_text <- paste0(
    "Panel ", panel_letters[[model_index]], ": ", specification
  )
  panel_lines <- strwrap(panel_text, width = 38)
  panel_heading <- paste0(
    "\\makecell[l]{",
    paste0(
      "\\textit{\\textbf{", latex_escape(panel_lines), "}}",
      collapse = "\\\\"
    ),
    "}"
  )
  latex_rows <- c(
    latex_rows,
    "\\hline",
    paste0(
      "\\multicolumn{4}{@{}l}{\\makebox[0pt][l]{",
      panel_heading,
      "}}\\\\"
    ),
    "\\hline"
  )
  for (treatment in required_treatments) {
    cells <- character(length(target_distances))
    for (distance_index in seq_along(target_distances)) {
      row <- selected$Specification == specification &
        selected$Treatment == treatment &
        selected$target_distance == target_distances[[distance_index]]
      if (sum(row) != 1L) {
        stop("Could not uniquely locate a multiplier table cell.", call. = FALSE)
      }
      cells[[distance_index]] <- sprintf(
        "\\makecell[c]{%s\\\\(%s, %s)}",
        format_number(selected$median[row]),
        format_number(selected$lower[row]),
        format_number(selected$upper[row])
      )
    }
    latex_rows <- c(
      latex_rows,
      paste0(
        "\\hspace{1em}", latex_escape(treatment), " & ",
        paste(cells, collapse = " & "), "\\\\"
      ),
      "\\addlinespace"
    )
  }
}
latex <- c(
  "\\centering",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "\\multicolumn{4}{c}{Structural Social Multiplier Robustness} \\\\",
  "\\addlinespace",
  "Treatment & $M(500\\mathrm{m})$ & $M(1500\\mathrm{m})$ & $M(2500\\mathrm{m})$ \\\\",
  "\\midrule",
  latex_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}"
)
dir.create(dirname(tex_output), recursive = TRUE, showWarnings = FALSE)
writeLines(latex, tex_output)

message("Wrote social-multiplier robustness summary: ", csv_output)
message("Wrote stripped Overleaf table: ", tex_output)
