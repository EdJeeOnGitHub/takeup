library(tidyverse)
library(fixest)

params <- lst(
  table_output_path = "presentations/rf-tables/main-specs"
)

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source("analysis_util.R")
source("dist_structural_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("scratch", "reduced-form-setup.R"))

dir.create(params$table_output_path, showWarnings = FALSE, recursive = TRUE)

treat_terms <- c(
  "assigned_treatment::ink",
  "assigned_treatment::calendar",
  "assigned_treatment::bracelet"
)

dist_treat_terms <- c(
  "assigned_treatment::ink:assigned_dist_group::close",
  "assigned_treatment::ink:assigned_dist_group::far",
  "assigned_treatment::calendar:assigned_dist_group::close",
  "assigned_treatment::calendar:assigned_dist_group::far",
  "assigned_treatment::bracelet:assigned_dist_group::close",
  "assigned_treatment::bracelet:assigned_dist_group::far"
)

pval_fmt <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "$<$0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

mean_fmt <- function(x) sprintf("%.3f", x)

wald_p <- function(model, terms) {
  terms <- intersect(terms, names(coef(model)))
  if (length(terms) == 0) return("")
  pval_fmt(wald(model, keep = paste(terms, collapse = "|"))$p)
}

tex_postprocessing <- function(tex) {
  tex <- tex %>%
    str_remove("\\\\begin\\{table\\}\\[htbp\\]") %>%
    str_remove("\\\\end\\{table\\}") %>%
    str_replace("(?m)^(\\s*)Female", "\\1\\\\addlinespace[0.4em]\n\\1Female")

  lines <- str_split(tex, "\n", simplify = FALSE)[[1]]
  female_row <- which(str_detect(lines, "^\\s*Female\\b"))[1]
  distance_row <- which(str_detect(lines, "^\\s*Expected distance to PoT\\b"))[1]
  control_mean_row <- which(str_detect(lines, "^\\s*Control mean\\b"))[1]

  if (any(is.na(c(female_row, distance_row, control_mean_row)))) {
    return(tex)
  }

  cov_start <- if (female_row > 1 && str_detect(lines[female_row - 1], "\\\\addlinespace")) {
    female_row - 1
  } else {
    female_row
  }
  cov_end <- min(distance_row + 1, length(lines))

  if (cov_start > control_mean_row) {
    return(tex)
  }

  cov_block <- lines[cov_start:cov_end]
  remaining <- lines[-c(cov_start:cov_end)]
  insert_before <- which(str_detect(remaining, "^\\s*Control mean\\b"))[1]

  paste(c(
    remaining[seq_len(insert_before - 1)],
    cov_block,
    remaining[insert_before:length(remaining)]
  ), collapse = "\n")
}

move_covariates_to_bottom <- function(path) {
  lines <- readLines(path, warn = FALSE)
  female_row <- which(str_detect(lines, "^\\s*Female\\b"))[1]
  distance_row <- which(str_detect(lines, "^\\s*Expected distance to PoT\\b"))[1]
  control_mean_row <- which(str_detect(lines, "^\\s*Control mean\\b"))[1]

  if (any(is.na(c(female_row, distance_row, control_mean_row)))) {
    return(invisible(FALSE))
  }

  cov_start <- if (female_row > 1 && str_detect(lines[female_row - 1], "\\\\addlinespace")) {
    female_row - 1
  } else {
    female_row
  }
  cov_end <- min(distance_row + 1, length(lines))

  if (cov_start > control_mean_row) {
    return(invisible(FALSE))
  }

  cov_block <- lines[cov_start:cov_end]
  remaining <- lines[-c(cov_start:cov_end)]
  insert_before <- which(str_detect(remaining, "^\\s*Control mean\\b"))[1]

  lines_out <- c(
    remaining[seq_len(insert_before - 1)],
    cov_block,
    remaining[insert_before:length(remaining)]
  )
  writeLines(lines_out, path)
  invisible(TRUE)
}

finalize_overleaf_input_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!str_detect(lines, "^\\s*\\\\begingroup\\s*$")]
  lines <- lines[!str_detect(lines, "^\\s*\\\\par\\\\endgroup\\s*$")]
  lines <- lines[!str_detect(lines, "^\\s*\\\\caption\\{")]
  lines <- lines[!str_detect(lines, "Clustered \\(cluster\\.id\\) standard-errors")]
  lines <- lines[!str_detect(lines, "Signif\\. Codes:")]
  lines <- lines[!str_detect(lines, "^\\s*\\\\emph\\{Variables\\}")]
  lines <- lines[!str_detect(lines, "^\\s*\\\\emph\\{Fit statistics\\}")]
  lines <- str_replace(
    lines,
    "^\\s*Model:\\s*&",
    "   Dependent variable: Missing from monitored take-up sample       &"
  )

  par_row <- which(str_detect(lines, "^\\\\par \\\\raggedright"))[1]
  if (!is.na(par_row)) {
    lines <- lines[seq_len(par_row - 1)]
  }

  tabular_row <- which(str_detect(lines, "^\\\\begin\\{tabular\\}"))[1]
  if (!is.na(tabular_row)) {
    lines <- append(
      lines,
      "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{",
      after = tabular_row - 1
    )
  }

  end_row <- which(str_detect(lines, "^\\\\end\\{tabular\\}"))[1]
  if (!is.na(end_row)) {
    lines[end_row] <- "\\end{tabular}}"
  }

  writeLines(lines, path)
  invisible(TRUE)
}

cluster_controls <- cov_analysis_data %>%
  select(cluster.id.x, mu_d, standard_clust_expected_dist, standard_cluster.dist.to.pot) %>%
  mutate(cluster.id = as.numeric(cluster.id.x)) %>%
  distinct(cluster.id, .keep_all = TRUE)

dropped_rows <- census_data %>%
  filter(!is.na(cluster.id)) %>%
  filter(sms.treatment == "sms.control") %>%
  filter(monitored) %>%
  filter(have_phone == "No" | have_phone == "Don't know number") %>%
  mutate(dropped_cto = monitored & !true.monitored) %>%
  filter(dropped_cto) %>%
  transmute(
    KEY.individ,
    cluster.id = as.numeric(cluster.id),
    dropped_cto,
    female = case_when(
      gender == "female" ~ TRUE,
      gender == "male" ~ FALSE,
      gender == 2 ~ TRUE,
      gender == 1 ~ FALSE,
      TRUE ~ NA
    ),
    age.census
  ) %>%
  left_join(
    analysis_data %>%
      transmute(
        cluster.id = as.numeric(levels(cluster.id)[cluster.id]),
        county,
        assigned_treatment,
        assigned_dist_group,
        cluster_id,
        cluster_id_rank
      ) %>%
      distinct(),
    by = "cluster.id"
  ) %>%
  left_join(cluster_controls, by = "cluster.id")

observed_rows <- analysis_data %>%
  transmute(
    KEY.individ,
    cluster.id = as.numeric(levels(cluster.id)[cluster.id]),
    dropped_cto = FALSE,
    female,
    age.census,
    county,
    assigned_treatment,
    assigned_dist_group,
    cluster_id,
    cluster_id_rank,
    standard_cluster.dist.to.pot
  ) %>%
  left_join(cluster_controls, by = "cluster.id", suffix = c("", ".cluster")) %>%
  mutate(
    standard_cluster.dist.to.pot = coalesce(
      standard_cluster.dist.to.pot,
      standard_cluster.dist.to.pot.cluster
    )
  ) %>%
  select(-ends_with(".cluster"))

main_analysis_attrition_df <- bind_rows(dropped_rows, observed_rows) %>%
  mutate(
    dropped_cto = as.numeric(dropped_cto),
    assigned_treatment = factor(as.character(assigned_treatment), levels = c("control", "ink", "calendar", "bracelet")),
    assigned_dist_group = factor(as.character(assigned_dist_group), levels = c("close", "far")),
    county = factor(county),
    cluster.id = factor(cluster.id)
  )

stopifnot(nrow(main_analysis_attrition_df) == 10794)

fit_treat <- feols(
  dropped_cto ~ i(assigned_treatment, ref = "control") +
    female + age.census + mu_d | county,
  data = main_analysis_attrition_df,
  cluster = ~cluster.id
)

fit_dist_controls <- feols(
  dropped_cto ~ assigned_dist_group +
    i(assigned_treatment, assigned_dist_group, ref = "control") +
    female + age.census + mu_d | county,
  data = main_analysis_attrition_df,
  cluster = ~cluster.id
)

control_mean <- main_analysis_attrition_df %>%
  filter(assigned_treatment == "control") %>%
  summarise(mean = mean(dropped_cto)) %>%
  pull(mean) %>%
  mean_fmt()

setFixest_dict(c(
  dropped_cto = "Dropped from monitored sample",
  "assigned_treatment::ink" = "Ink",
  "assigned_treatment::calendar" = "Calendar",
  "assigned_treatment::bracelet" = "Bracelet",
  assigned_dist_groupfar = "Far",
  "assigned_treatment::ink:assigned_dist_group::close" = "Ink $\\times$ Close",
  "assigned_treatment::ink:assigned_dist_group::far" = "Ink $\\times$ Far",
  "assigned_treatment::calendar:assigned_dist_group::close" = "Calendar $\\times$ Close",
  "assigned_treatment::calendar:assigned_dist_group::far" = "Calendar $\\times$ Far",
  "assigned_treatment::bracelet:assigned_dist_group::close" = "Bracelet $\\times$ Close",
  "assigned_treatment::bracelet:assigned_dist_group::far" = "Bracelet $\\times$ Far",
  femaleTRUE = "Female",
  age.census = "Age",
  mu_d = "Expected distance to PoT"
))

table_file <- file.path(params$table_output_path, "paper_main_analysis_sample_attrition_tbl.tex")

etable(
  fit_treat,
  fit_dist_controls,
  headers = c("Treatment", "Treatment $\\times$ Distance"),
  depvar = FALSE,
  fitstat = c("n", "r2"),
  se.below = TRUE,
  tex = TRUE,
  file = table_file,
  replace = TRUE,
  postprocess.tex = tex_postprocessing,
  digits = 3,
  digits.stats = 3,
  drop.section = "fixef",
  order = c(
    "assigned_treatment::ink",
    "assigned_treatment::calendar",
    "assigned_treatment::bracelet",
    "assigned_dist_groupfar",
    "assigned_treatment::ink:assigned_dist_group::close",
    "assigned_treatment::ink:assigned_dist_group::far",
    "assigned_treatment::calendar:assigned_dist_group::close",
    "assigned_treatment::calendar:assigned_dist_group::far",
    "assigned_treatment::bracelet:assigned_dist_group::close",
    "assigned_treatment::bracelet:assigned_dist_group::far",
    "femaleTRUE",
    "age.census",
    "mu_d"
  ),
  style.df = style.df(depvar.title = "", fixef.title = "", var.title = "", stats.title = ""),
  extralines = list(
    "Control mean" = rep(control_mean, 2),
    "County FEs" = rep("Yes", 2),
    "RF controls" = rep("Yes", 2),
    "$H_0$: all displayed treatment coefficients = 0, $p$-value" = c(
      wald_p(fit_treat, treat_terms),
      wald_p(fit_dist_controls, dist_treat_terms)
    )
  )
)

move_covariates_to_bottom(table_file)
finalize_overleaf_input_table(table_file)
