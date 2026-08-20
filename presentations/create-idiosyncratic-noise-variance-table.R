#!/usr/bin/env Rscript

read_stan_csv_column <- function(files, colname) {
  unlist(lapply(files, function(file) {
    if (!file.exists(file)) {
      stop("Missing Stan CSV: ", file)
    }

    cmd <- sprintf(
      "awk -F, 'BEGIN{idx=0} /^#/ {next} idx==0 {for(i=1;i<=NF;i++) if($i==\"%s\") idx=i; if(idx==0) exit 2; next} {print $idx}' %s",
      colname,
      shQuote(file)
    )
    con <- pipe(cmd)
    on.exit(close(con), add = TRUE)
    scan(con, quiet = TRUE)
  }), use.names = FALSE)
}

summarise_sigma_var <- function(u_sd) {
  sigma_var <- u_sd^2
  qs <- quantile(sigma_var, c(0.025, 0.5, 0.975), na.rm = TRUE)
  c(
    mean = mean(sigma_var, na.rm = TRUE),
    median = qs[[2]],
    lo = qs[[1]],
    hi = qs[[3]]
  )
}

format_num <- function(x) sprintf("%.3f", x)
format_pct <- function(x) sprintf("%+.1f\\%%", 100 * x)

main_files <- file.path(
  "data/stan_analysis_data",
  "dist_fit95_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB-1.csv"
)

main_stats <- summarise_sigma_var(read_stan_csv_column(main_files, "u_sd.1"))

# Computed by streaming u_sd.1 from the four fit105
# STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP Stan CSV files on
# Midway. The raw CSVs are not stored locally in this checkout.
indiv_stats <- c(
  mean = 0.38030,
  median = 0.16108,
  lo = 0.019920,
  hi = 2.11670
)

mean_diff <- indiv_stats[["mean"]] / main_stats[["mean"]] - 1
median_diff <- indiv_stats[["median"]] / main_stats[["median"]] - 1

table_tex <- sprintf(
  paste0(
    "\n",
    "\\centering\n\n",
    "\\centering\n",
    "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{\n",
    "\\begin{tabular}[t]{lcccc}\n",
    "\\toprule\n",
    "Specification & Mean & Median & 95\\%% credible interval & Difference vs. main\\\\\n",
    "\\midrule\n",
    "Main specification & %s & %s & (%s, %s) & --\\\\\n",
    "Individual distance, individual fixed point & %s & %s & (%s, %s) & \\makecell[c]{%s mean\\\\%s median}\\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}}\n"
  ),
  format_num(main_stats[["mean"]]),
  format_num(main_stats[["median"]]),
  format_num(main_stats[["lo"]]),
  format_num(main_stats[["hi"]]),
  format_num(indiv_stats[["mean"]]),
  format_num(indiv_stats[["median"]]),
  format_num(indiv_stats[["lo"]]),
  format_num(indiv_stats[["hi"]]),
  format_pct(mean_diff),
  format_pct(median_diff)
)

output_file <- "presentations/tables/fit105/idiosyncratic-noise-variance-comparison.tex"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
writeLines(table_tex, output_file)
