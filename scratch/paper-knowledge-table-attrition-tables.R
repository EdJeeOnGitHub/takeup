library(tidyverse)
library(fixest)

set.seed(20260710)

params <- lst(
  table_output_path = "presentations/rf-tables/main-specs",
  output_path = "temp-data",
  n_sims = 1000
)

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source("analysis_util.R")
source("dist_structural_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("scratch", "reduced-form-setup.R"))

dir.create(params$table_output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(params$output_path, showWarnings = FALSE, recursive = TRUE)

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

cov_terms <- c("femaleTRUE", "age.census", "mu_d")

pval_fmt <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "$<$0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

fmt <- function(x) sprintf("%.3f", x)

stars <- function(p) case_when(
  is.na(p) ~ "",
  p < 0.01 ~ "$^{***}$",
  p < 0.05 ~ "$^{**}$",
  p < 0.1 ~ "$^{*}$",
  TRUE ~ ""
)

coef_cell <- function(estimate, se, pval) {
  if (is.na(estimate)) return("")
  paste0(fmt(estimate), stars(pval), "\\\\{(", fmt(se), ")}")
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
    "   Dependent variable: Missing from knowledge table                &"
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

wald_p <- function(model, terms) {
  terms <- intersect(terms, names(coef(model)))
  if (length(terms) == 0) return("")
  pval_fmt(wald(model, keep = paste(terms, collapse = "|"))$p)
}

setFixest_dict(c(
  attrit_know_table = "Missing knowledge table",
  "assigned_treatment::ink" = "Ink",
  "assigned_treatment::calendar" = "Calendar",
  "assigned_treatment::bracelet" = "Bracelet",
  "assigned_treatment::ink:assigned_dist_group::close" = "Ink $\\times$ Close",
  "assigned_treatment::ink:assigned_dist_group::far" = "Ink $\\times$ Far",
  "assigned_treatment::calendar:assigned_dist_group::close" = "Calendar $\\times$ Close",
  "assigned_treatment::calendar:assigned_dist_group::far" = "Calendar $\\times$ Far",
  "assigned_treatment::bracelet:assigned_dist_group::close" = "Bracelet $\\times$ Close",
  "assigned_treatment::bracelet:assigned_dist_group::far" = "Bracelet $\\times$ Far",
  assigned_dist_groupfar = "Far",
  femaleTRUE = "Female",
  age.census = "Age",
  mu_d = "Expected distance to PoT"
))

full_sample_attrition_df <- all_endline_data %>%
  left_join(
    all_data %>%
      select(KEY.individ, all_of(l_cov_vars)),
    by = "KEY.individ"
  ) %>%
  left_join(
    cov_analysis_data %>%
      select(cluster.id.x, mu_d, standard_clust_expected_dist, standard_cluster.dist.to.pot) %>%
      mutate(cluster.id = as.numeric(cluster.id.x)) %>%
      unique(),
    by = "cluster.id"
  ) %>%
  mutate(
    attrit_know_table = !in_know_table,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group
  )

full_fit_treat <- feols(
  attrit_know_table ~ i(assigned_treatment, ref = "control") +
    .[l_cov_vars] + mu_d | county,
  data = full_sample_attrition_df,
  cluster = ~cluster.id
)

full_fit_dist <- feols(
  attrit_know_table ~ assigned_dist_group +
    i(assigned_treatment, assigned_dist_group, ref = "control") +
    .[l_cov_vars] + mu_d | county,
  data = full_sample_attrition_df,
  cluster = ~cluster.id
)

full_control_mean <- full_sample_attrition_df %>%
  filter(assigned_treatment == "control") %>%
  summarise(mean = mean(attrit_know_table)) %>%
  pull(mean) %>%
  fmt()

full_table_file <- file.path(params$table_output_path, "paper_full_knowledge_table_attrition_tbl.tex")

etable(
  full_fit_treat,
  full_fit_dist,
  headers = c("Treatment", "Treatment $\\times$ Distance"),
  depvar = FALSE,
  fitstat = c("n", "r2"),
  se.below = TRUE,
  tex = TRUE,
  file = full_table_file,
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
    "Control mean" = rep(full_control_mean, 2),
    "County FEs" = rep("Yes", 2),
    "RF controls" = rep("Yes", 2),
    "$H_0$: all displayed treatment coefficients = 0, $p$-value" = c(
      wald_p(full_fit_treat, treat_terms),
      wald_p(full_fit_dist, dist_treat_terms)
    )
  )
)

move_covariates_to_bottom(full_table_file)
finalize_overleaf_input_table(full_table_file)

summ_know_A_df <- summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A")) %>%
  select(KEY.individ, know.table.type, obs_know_person)

base_df <- endline_data %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  mutate(
    observed_table_a = !is.na(know.table.type) & obs_know_person > 0,
    missing_know_table = !in_know_table,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group
  )

observed_a_df <- base_df %>%
  filter(observed_table_a) %>%
  mutate(attrit_know_table = FALSE)

missing_df <- base_df %>%
  filter(missing_know_table) %>%
  mutate(attrit_know_table = TRUE)

stopifnot(nrow(observed_a_df) == 1141)
stopifnot(nrow(missing_df) == 252)

model_terms <- list(
  treat_controls = c(treat_terms, cov_terms),
  dist_controls = c("assigned_dist_groupfar", dist_treat_terms, cov_terms)
)

model_test_terms <- list(
  treat_controls = treat_terms,
  dist_controls = dist_treat_terms
)

fit_sim_models <- function(sim_df) {
  list(
    treat_controls = feols(
      attrit_know_table ~ i(assigned_treatment, ref = "control") +
        .[l_cov_vars] + mu_d | county,
      data = sim_df,
      cluster = ~cluster.id
    ),
    dist_controls = feols(
      attrit_know_table ~ assigned_dist_group +
        i(assigned_treatment, assigned_dist_group, ref = "control") +
        .[l_cov_vars] + mu_d | county,
      data = sim_df,
      cluster = ~cluster.id
    )
  )
}

extract_sim_model <- function(model, sim, model_name) {
  terms <- model_terms[[model_name]]
  coefs <- coef(model)
  vc <- vcov(model)
  present <- intersect(terms, names(coefs))

  coef_rows <- map_dfr(present, function(term) {
    tibble(
      sim = sim,
      model = model_name,
      term = term,
      estimate = unname(coefs[term]),
      se = sqrt(vc[term, term])
    )
  })

  vc_rows <- expand_grid(term1 = present, term2 = present) %>%
    mutate(
      sim = sim,
      model = model_name,
      vcov = map2_dbl(term1, term2, ~ vc[.x, .y])
    )

  list(coefs = coef_rows, vcov = vc_rows)
}

fit_one_sim <- function(sim) {
  imputed_missing_a <- missing_df %>%
    slice_sample(n = nrow(missing_df) / 2)

  sim_df <- bind_rows(observed_a_df, imputed_missing_a)
  models <- fit_sim_models(sim_df)

  extracted <- imap(models, ~ extract_sim_model(.x, sim, .y))

  control_mean <- sim_df %>%
    filter(assigned_treatment == "control") %>%
    summarise(control_mean = mean(attrit_know_table), .groups = "drop") %>%
    mutate(sim = sim)

  list(
    coefs = map_dfr(extracted, "coefs"),
    vcov = map_dfr(extracted, "vcov"),
    control = control_mean
  )
}

sim_coef_path <- file.path(params$output_path, "paper-simulated-table-A-attrition-coefs.csv")
sim_vcov_path <- file.path(params$output_path, "paper-simulated-table-A-attrition-vcov.csv")
sim_control_path <- file.path(params$output_path, "paper-simulated-table-A-attrition-control-means.csv")

if (file.exists(sim_coef_path) && file.exists(sim_vcov_path) && file.exists(sim_control_path)) {
  sim_coefs <- read_csv(sim_coef_path, show_col_types = FALSE)
  sim_vcov <- read_csv(sim_vcov_path, show_col_types = FALSE)
  sim_control <- read_csv(sim_control_path, show_col_types = FALSE)
} else {
  sim_results <- map(seq_len(params$n_sims), fit_one_sim, .progress = TRUE)
  sim_coefs <- map_dfr(sim_results, "coefs")
  sim_vcov <- map_dfr(sim_results, "vcov")
  sim_control <- map_dfr(sim_results, "control")

  write_csv(sim_coefs, sim_coef_path)
  write_csv(sim_vcov, sim_vcov_path)
  write_csv(sim_control, sim_control_path)
}

combine_model <- function(model_name) {
  model_coefs <- sim_coefs %>%
    filter(model == model_name)
  model_vcov <- sim_vcov %>%
    filter(model == model_name)

  terms <- model_terms[[model_name]]
  present <- intersect(terms, unique(model_coefs$term))
  m <- n_distinct(model_coefs$sim)

  q_mat <- model_coefs %>%
    filter(term %in% present) %>%
    select(sim, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    arrange(sim)

  q <- as.matrix(q_mat[, present, drop = FALSE])
  q_bar <- colMeans(q, na.rm = TRUE)

  u_bar <- matrix(0, nrow = length(present), ncol = length(present), dimnames = list(present, present))
  for (s in q_mat$sim) {
    vc_s <- model_vcov %>%
      filter(sim == s, term1 %in% present, term2 %in% present) %>%
      select(term1, term2, vcov) %>%
      pivot_wider(names_from = term2, values_from = vcov)
    vc_mat <- as.matrix(vc_s[match(present, vc_s$term1), present, drop = FALSE])
    u_bar <- u_bar + vc_mat / m
  }

  b <- stats::cov(q)
  total_v <- u_bar + (1 + 1 / m) * b

  tibble(
    model = model_name,
    term = present,
    estimate = unname(q_bar[present]),
    se = sqrt(diag(total_v)[present]),
    pval = 2 * pnorm(-abs(unname(q_bar[present]) / sqrt(diag(total_v)[present])))
  ) %>%
    bind_rows(
      tibble(
        model = model_name,
        term = "joint_treatment_p",
        estimate = NA_real_,
        se = NA_real_,
        pval = {
          test_terms <- intersect(model_test_terms[[model_name]], present)
          q_test <- q_bar[test_terms]
          v_test <- total_v[test_terms, test_terms, drop = FALSE]
          stat <- drop(t(q_test) %*% solve(v_test) %*% q_test)
          pchisq(stat, df = length(test_terms), lower.tail = FALSE)
        }
      )
    )
}

sim_combined <- map_dfr(names(model_terms), combine_model)
write_csv(sim_combined, file.path(params$output_path, "paper-simulated-table-A-attrition-combined.csv"))

sim_control_mean <- mean(sim_control$control_mean)
sim_control_se <- sqrt(var(sim_control$control_mean) * (1 + 1 / nrow(sim_control)))

sim_cell <- function(model, term) {
  row <- sim_combined %>% filter(model == !!model, term == !!term)
  if (nrow(row) == 0) return("")
  coef_cell(row$estimate, row$se, row$pval)
}

sim_pval_cell <- function(model) {
  row <- sim_combined %>% filter(model == !!model, term == "joint_treatment_p")
  if (nrow(row) == 0) return("")
  pval_fmt(row$pval)
}

rows <- tribble(
  ~label, ~term, ~models,
  "Ink", "assigned_treatment::ink", list(c("treat_controls")),
  "Calendar", "assigned_treatment::calendar", list(c("treat_controls")),
  "Bracelet", "assigned_treatment::bracelet", list(c("treat_controls")),
  "Far", "assigned_dist_groupfar", list(c("dist_controls")),
  "Ink $\\times$ Close", "assigned_treatment::ink:assigned_dist_group::close", list(c("dist_controls")),
  "Ink $\\times$ Far", "assigned_treatment::ink:assigned_dist_group::far", list(c("dist_controls")),
  "Calendar $\\times$ Close", "assigned_treatment::calendar:assigned_dist_group::close", list(c("dist_controls")),
  "Calendar $\\times$ Far", "assigned_treatment::calendar:assigned_dist_group::far", list(c("dist_controls")),
  "Bracelet $\\times$ Close", "assigned_treatment::bracelet:assigned_dist_group::close", list(c("dist_controls")),
  "Bracelet $\\times$ Far", "assigned_treatment::bracelet:assigned_dist_group::far", list(c("dist_controls")),
  "Female", "femaleTRUE", list(c("treat_controls", "dist_controls")),
  "Age", "age.census", list(c("treat_controls", "dist_controls")),
  "Expected distance to PoT", "mu_d", list(c("treat_controls", "dist_controls"))
)

sim_coef_rows <- rows %>%
  rowwise() %>%
  mutate(
    c1 = if ("treat_controls" %in% unlist(models)) sim_cell("treat_controls", term) else "",
    c2 = if ("dist_controls" %in% unlist(models)) sim_cell("dist_controls", term) else "",
    row = paste0(label, " & \\makecell[c]{", c1, "} & \\makecell[c]{", c2,
                 "} \\\\")
  ) %>%
  ungroup() %>%
  pull(row)

sim_coef_rows <- str_replace(
  sim_coef_rows,
  "^Female",
  "\\\\addlinespace[0.4em]\nFemale"
)

sim_tex_table <- c(
  "\\centering",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{",
  "\\begin{tabular}[t]{lcc}",
  "\\toprule",
  " & Treatment & Treatment $\\times$ Distance \\\\",
  "Dependent variable: Missing from observability module & (1) & (2) \\\\",
  "\\midrule",
  sim_coef_rows,
  "\\midrule",
  paste0("Control mean & ", fmt(sim_control_mean), " & ", fmt(sim_control_mean), " \\\\"),
  "County FEs & Yes & Yes \\\\",
  "RF controls & Yes & Yes \\\\",
  paste0("$H_0$: all displayed treatment coefficients = 0, $p$-value & ",
         sim_pval_cell("treat_controls"), " & ", sim_pval_cell("dist_controls"), " \\\\"),
  "Observations & 1,267 & 1,267 \\\\",
  "Simulations & 1,000 & 1,000 \\\\",
  "\\bottomrule",
  "\\end{tabular}}"
)

writeLines(sim_tex_table, file.path(params$table_output_path, "paper_simulated_table_A_attrition_tbl.tex"))

preview_tex <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{booktabs}",
  "\\usepackage{graphicx}",
  "\\newcommand{\\makecell}[2][]{\\begin{tabular}{@{}c@{}}#2\\end{tabular}}",
  "\\begin{document}",
  "\\begin{table}[htbp]",
  "\\scriptsize",
  "\\input{../presentations/rf-tables/main-specs/paper_full_knowledge_table_attrition_tbl.tex}",
  "\\end{table}",
  "\\clearpage",
  "\\begin{table}[htbp]",
  "\\scriptsize",
  "\\input{../presentations/rf-tables/main-specs/paper_simulated_table_A_attrition_tbl.tex}",
  "\\end{table}",
  "\\end{document}"
)

writeLines(preview_tex, "scratch/paper-knowledge-table-attrition-preview.tex")
