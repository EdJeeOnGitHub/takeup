library(tidyverse)
library(broom)
library(data.table)
library(kableExtra)
library(knitr)
library(fixest)

script_options <- docopt::docopt(
  "Usage:
  reduced-form-bootstrap.R [options]

Options:
  --plots              Run NP fit plots section
  --externality        Run externality knowledge section
  --takeup             Run takeup regressions + levels section
  --itt                Run original-distance-assignment take-up ITT robustness
  --beliefs            Run beliefs (FOB/SOB) regressions + levels section
  --endline            Run endline/incentive/preference/travel section
  --sms                Run SMS heterogeneity section
  --heterogeneity      Run heterogeneity + WTP section
  --table-output-path=<path>  Table output path [default: presentations/rf-tables/main-specs]
  --stat=<stat>        Statistic to show [default: std.error]
",
  args = if (interactive()) "
   --endline 
  " else commandArgs(trailingOnly = TRUE)
)
# TODO: Fix --sms and --heterogeneity

run_all <- !any(script_options$plots, script_options$externality, script_options$takeup,
                script_options$itt,
                script_options$beliefs, script_options$endline, script_options$sms,
                script_options$heterogeneity)

tex_postprocessing = function(tex) {
    tex %>%
        str_remove("\\\\begin\\{table\\}\\[htbp\\]") %>%
        str_remove("\\\\end\\{table\\}") %>%
        str_replace(
          .,
          "Covariate",
          "\\\\midrule\n   Covariate"
        )
}

drop_fixest_footer = function(tex) {
  tex %>%
    str_remove_all("(?m)^\\s*\\\\multicolumn\\{\\d+\\}\\{l\\}\\{\\\\emph\\{Clustered.*\\}\\}\\\\\\\\\\s*\n") %>%
    str_remove_all("(?m)^\\s*\\\\multicolumn\\{\\d+\\}\\{l\\}\\{\\\\emph\\{Signif\\. Codes:.*\\}\\}\\\\\\\\\\s*\n")
}

compact_fixest_postprocessing = function(tex) {
  tex %>%
    tex_postprocessing() %>%
    drop_fixest_footer()
}
# Prefix all messages/output with current section and report errors
current_section <- "(setup)"
section_message <- function(...) message(sprintf("[%s] %s", current_section, paste0(...)))
run_section <- function(name, expr) {
  current_section <<- name
  section_message("Starting")
  tryCatch(
    {
      result <- force(expr)
      section_message("Done")
      invisible(result)
    },
    error = function(e) {
      section_message("ERROR: ", conditionMessage(e))
      stop(sprintf("[%s] %s", name, conditionMessage(e)), call. = FALSE)
    }
  )
}

params <- lst(
  table_output_path = script_options$table_output_path,
  show_probs = FALSE,
  width = 0.95,
  cache = FALSE,
  fit = FALSE,
  stat = script_options$stat
)

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path("dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("scratch", "reduced-form-setup.R"))
# From running:
# pdslasso dewormed_num dpf ($cov_vars i.county_fac mu_d), cluster(clusteridx) pnotpen(i.county_fac)
# where mu_d is the expected distance to the cluster
# get the same result using actual distance.
# using treatment dummies and putting them in the amelioration set includes everything 
# which seems wrong
l_cov_vars = c(
  "female",
  "age.census"
)

assert_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

clean_cluster_id <- function(x, label = "cluster ID") {
  x_chr <- str_trim(as.character(x))
  bad <- is.na(x_chr) | x_chr == "" | !str_detect(x_chr, "^[0-9]+$")
  if (any(bad)) {
    bad_examples <- unique(x_chr[bad])
    bad_examples <- bad_examples[seq_len(min(length(bad_examples), 5))]
    stop(
      sprintf(
        "%s contains non-integer or missing values: %s",
        label,
        paste(bad_examples, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  as.integer(x_chr)
}

assert_no_cluster_id_collapse <- function(data, raw_col, clean_col, source_name) {
  id_map <- data %>%
    distinct(raw_id = .data[[raw_col]], clean_id = .data[[clean_col]])

  assert_true(
    !any(is.na(id_map$raw_id)) && !any(is.na(id_map$clean_id)),
    paste(source_name, "has missing cluster IDs after cleaning.")
  )
  assert_true(
    id_map %>% count(clean_id) %>% filter(n > 1) %>% nrow() == 0,
    paste(source_name, "has multiple raw cluster IDs mapping to one cleaned cluster ID.")
  )
}



#### Original Distance Assignment ITT ------------------------------------------

if (script_options$itt) run_section("Original Distance Assignment ITT", {

original_distance_df <- read_rds("data/rct_targetable_schools_2.0.rds") %>%
  as_tibble() %>%
  transmute(
    original_cluster_id_raw = as.character(pot.cluster.id),
    original_cluster_id = clean_cluster_id(pot.cluster.id, "rct_targetable_schools_2.0$pot.cluster.id"),
    original_dist_group = factor(assigned.dist.cat, levels = c("close", "far"))
  ) %>%
  distinct()

assert_no_cluster_id_collapse(
  original_distance_df,
  "original_cluster_id_raw",
  "original_cluster_id",
  "Original distance assignment data"
)
assert_true(
  nrow(original_distance_df) == n_distinct(original_distance_df$original_cluster_id),
  "Original distance assignments are not unique by original cluster ID."
)
assert_true(!any(is.na(original_distance_df$original_dist_group)), "Original distance assignments contain labels other than close/far.")

target_village_distance_df <- read_rds("data/rct_target_villages_2.0.rds") %>%
  as_tibble() %>%
  transmute(
    target_village_cluster_id_raw = as.character(cluster.id),
    original_cluster_id = clean_cluster_id(cluster.id, "rct_target_villages_2.0$cluster.id"),
    target_village_dist_group = factor(village.dist.cat, levels = c("close", "far"))
  ) %>%
  distinct()

assert_no_cluster_id_collapse(
  target_village_distance_df,
  "target_village_cluster_id_raw",
  "original_cluster_id",
  "Target-village distance data"
)
assert_true(
  nrow(target_village_distance_df) == n_distinct(target_village_distance_df$original_cluster_id),
  "Target-village distance assignments are not unique by original cluster ID."
)

original_distance_itt_data <- cov_analysis_data %>%
  mutate(
    analysis_cluster_id_raw = as.character(cluster.id.x),
    original_cluster_id = clean_cluster_id(cluster.id.x, "analysis cluster.id.x"),
    final_dist_group = factor(dist.pot.group, levels = c("close", "far")),
    cluster_id_rank = dense_rank(cluster_id)
  ) %>%
  left_join(
    original_distance_df,
    by = "original_cluster_id",
    relationship = "many-to-one"
  ) %>%
  mutate(
    assigned_dist_group = factor(original_dist_group, levels = c("close", "far")),
    signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
    signal = factor(signal, levels = c("no signal", "signal"))
  )

assert_true(nrow(original_distance_itt_data) == nrow(cov_analysis_data), "Joining original distance assignments changed the row count.")
assert_true(!anyDuplicated(original_distance_itt_data$KEY.individ), "KEY.individ is duplicated after joining original distance assignments.")
assert_true(!any(is.na(original_distance_itt_data$assigned_dist_group)), "Missing original distance assignments for analysis clusters.")
assert_true(!any(is.na(original_distance_itt_data$mu_d)), "Missing expected-distance controls in original-distance ITT data.")
assert_true(
  original_distance_itt_data %>% distinct(original_cluster_id, cluster_id) %>% count(original_cluster_id) %>% filter(n > 1) %>% nrow() == 0,
  "An original cluster ID maps to multiple sequential cluster IDs."
)
assert_true(
  original_distance_itt_data %>% distinct(original_cluster_id, cluster_id) %>% count(cluster_id) %>% filter(n > 1) %>% nrow() == 0,
  "A sequential cluster ID maps to multiple original cluster IDs."
)
assert_no_cluster_id_collapse(
  original_distance_itt_data,
  "analysis_cluster_id_raw",
  "original_cluster_id",
  "Analysis data"
)

monitored_source <- read_rds("data/clean-data/monitored-nosms-takeup-data.rds") %>%
  transmute(
    KEY.individ,
    monitored_cluster_id_raw = as.character(cluster.id),
    source_cluster_id = clean_cluster_id(cluster.id, "monitored-nosms-takeup-data$cluster.id"),
    source_treatment = as.character(assigned.treatment),
    source_dist_group = as.character(dist.pot.group),
    source_dist_to_pot = dist.to.pot,
    source_cluster_dist_to_pot = cluster.dist.to.pot
  )

monitored_check <- original_distance_itt_data %>%
  select(
    KEY.individ,
    original_cluster_id,
    assigned.treatment,
    dist.pot.group,
    dist.to.pot,
    cluster.dist.to.pot
  ) %>%
  left_join(monitored_source, by = "KEY.individ", relationship = "one-to-one")

assert_true(!any(is.na(monitored_check$source_cluster_id)), "Some analysis rows did not match monitored source data.")
assert_true(all(monitored_check$original_cluster_id == monitored_check$source_cluster_id), "Original cluster IDs differ from monitored source cluster IDs.")
assert_true(all(as.character(monitored_check$assigned.treatment) == monitored_check$source_treatment), "Treatment assignments differ from monitored source data.")
assert_true(all(as.character(monitored_check$dist.pot.group) == monitored_check$source_dist_group), "Current distance groups differ from monitored source data.")
assert_true(max(abs(monitored_check$dist.to.pot - monitored_check$source_dist_to_pot), na.rm = TRUE) < 1e-8, "Individual PoT distances differ from monitored source data.")
assert_true(max(abs(monitored_check$cluster.dist.to.pot - monitored_check$source_cluster_dist_to_pot), na.rm = TRUE) < 1e-8, "Cluster PoT distances differ from monitored source data.")

processed_check <- original_distance_itt_data %>%
  distinct(original_cluster_id, assigned.treatment, dist.pot.group) %>%
  left_join(
    read_rds("data/takeup_processed_cluster_strat.rds") %>%
      transmute(
        processed_cluster_id_raw = as.character(cluster.id),
        original_cluster_id = clean_cluster_id(cluster.id, "takeup_processed_cluster_strat$cluster.id"),
        processed_treatment = as.character(assigned.treatment),
        processed_dist_group = as.character(dist.pot.group)
      ) %>%
      distinct(),
    by = "original_cluster_id",
    relationship = "one-to-one"
  )

assert_true(!any(is.na(processed_check$processed_treatment)), "Some clusters did not match processed randomization data.")
assert_true(all(as.character(processed_check$assigned.treatment) == processed_check$processed_treatment), "Treatment assignments differ from processed randomization data.")
assert_true(all(as.character(processed_check$dist.pot.group) == processed_check$processed_dist_group), "Current distance groups differ from processed randomization data.")

cluster_audit <- original_distance_itt_data %>%
  distinct(
    original_cluster_id,
    analysis_cluster_id_raw,
    cluster_id,
    assigned.treatment,
    final_dist_group,
    original_dist_group
  ) %>%
  left_join(
    original_distance_df %>% select(original_cluster_id, original_cluster_id_raw),
    by = "original_cluster_id",
    relationship = "many-to-one"
  ) %>%
  left_join(
    target_village_distance_df,
    by = "original_cluster_id",
    relationship = "many-to-one"
  ) %>%
  left_join(
    processed_check %>%
      distinct(original_cluster_id, processed_cluster_id_raw, processed_treatment, processed_dist_group),
    by = "original_cluster_id",
    relationship = "one-to-one"
  ) %>%
  mutate(
    target_matches_original = as.character(target_village_dist_group) == as.character(original_dist_group),
    processed_matches_analysis = processed_treatment == as.character(assigned.treatment) &
      processed_dist_group == as.character(final_dist_group),
    final_matches_original = as.character(final_dist_group) == as.character(original_dist_group),
    switch_direction = case_when(
      final_matches_original ~ "unchanged",
      TRUE ~ paste(original_dist_group, "to", final_dist_group)
    )
  ) %>%
  arrange(original_cluster_id)

write_csv(cluster_audit, "temp-data/original-distance-cluster-audit.csv")

assert_true(!any(is.na(cluster_audit$target_village_dist_group)), "Some analysis clusters did not match target-village data.")
assert_true(all(cluster_audit$target_matches_original), "Target-village distance assignments differ from original randomization assignments.")
assert_true(all(cluster_audit$processed_matches_analysis), "Processed randomization data differs from analysis data.")

merge_checks <- tibble(
  check = c(
    "analysis_rows",
    "analysis_clusters",
    "missing_original_distance_assignments",
    "missing_expected_distance_controls",
    "treatment_mismatches_monitored",
    "distance_group_mismatches_monitored",
    "max_abs_dist_to_pot_diff_monitored",
    "max_abs_cluster_dist_to_pot_diff_monitored",
    "treatment_mismatches_processed",
    "distance_group_mismatches_processed"
  ),
  value = c(
    nrow(original_distance_itt_data),
    n_distinct(original_distance_itt_data$original_cluster_id),
    sum(is.na(original_distance_itt_data$assigned_dist_group)),
    sum(is.na(original_distance_itt_data$mu_d)),
    sum(as.character(monitored_check$assigned.treatment) != monitored_check$source_treatment),
    sum(as.character(monitored_check$dist.pot.group) != monitored_check$source_dist_group),
    max(abs(monitored_check$dist.to.pot - monitored_check$source_dist_to_pot), na.rm = TRUE),
    max(abs(monitored_check$cluster.dist.to.pot - monitored_check$source_cluster_dist_to_pot), na.rm = TRUE),
    sum(as.character(processed_check$assigned.treatment) != processed_check$processed_treatment),
    sum(as.character(processed_check$dist.pot.group) != processed_check$processed_dist_group)
  )
)

write_csv(merge_checks, "temp-data/original-distance-merge-checks.csv")

original_distance_switches <- original_distance_itt_data %>%
  distinct(original_cluster_id, final_dist_group, original_dist_group) %>%
  count(original_dist_group, final_dist_group, name = "clusters") %>%
  arrange(original_dist_group, final_dist_group)

write_csv(original_distance_switches, "temp-data/original-distance-switch-counts.csv")

expected_original_distance_switches <- tribble(
  ~original_dist_group, ~final_dist_group, ~clusters,
  factor("close", levels = c("close", "far")), factor("close", levels = c("close", "far")), 64L,
  factor("close", levels = c("close", "far")), factor("far", levels = c("close", "far")), 10L,
  factor("far", levels = c("close", "far")), factor("close", levels = c("close", "far")), 16L,
  factor("far", levels = c("close", "far")), factor("far", levels = c("close", "far")), 54L
)

assert_true(
  setequal(
    original_distance_switches %>% mutate(across(c(original_dist_group, final_dist_group), as.character)),
    expected_original_distance_switches %>% mutate(across(c(original_dist_group, final_dist_group), as.character))
  ),
  "Original-distance switch counts differ from the documented 64/10/16/54 cluster split."
)

if (file.exists("docs/randomization-mismatches.csv")) {
  documented_switches <- read_csv("docs/randomization-mismatches.csv", show_col_types = FALSE) %>%
    transmute(
      original_cluster_id = clean_cluster_id(pot.cluster.id, "docs/randomization-mismatches.csv$pot.cluster.id"),
      documented_original_dist_group = as.character(clust_randomization_cat),
      documented_final_dist_group = as.character(processed_cat)
    ) %>%
    arrange(original_cluster_id)

  generated_switches <- cluster_audit %>%
    filter(!final_matches_original) %>%
    transmute(
      original_cluster_id,
      documented_original_dist_group = as.character(original_dist_group),
      documented_final_dist_group = as.character(final_dist_group)
    ) %>%
    arrange(original_cluster_id)

  assert_true(
    setequal(documented_switches, generated_switches),
    "Generated original-distance switches differ from docs/randomization-mismatches.csv."
  )
}

original_distance_discrete_regression <- function(data, weights) {
  feols(
    dewormed ~ 0 + assigned_treatment*assigned_dist_group + .[l_cov_vars] + mu_d | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

original_distance_itt_output <- wrapper_function(
  data = original_distance_itt_data,
  regression_spec = original_distance_discrete_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/original-distance-itt-tidy-tes.csv",
  table_name = "rf_original_distance_itt_tbl"
)

print(original_distance_switches)
original_distance_itt_output$tidy_summary %>%
  filter(assigned_dist_group == "far - close") %>%
  print(n = Inf)

summ_know_A_itt_df <- summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A"))

endline_know_A_itt_data <- endline_data %>%
  left_join(
    summ_know_A_itt_df,
    by = "KEY.individ"
  ) %>%
  filter(sms.treatment == "sms.control", obs_know_person > 0) %>%
  mutate(
    observability_cluster_id_raw = as.character(cluster.id),
    original_cluster_id = clean_cluster_id(cluster.id, "endline_data$cluster.id"),
    final_dist_group = factor(dist.pot.group, levels = c("close", "far")),
    cluster_id_rank = dense_rank(original_cluster_id),
    prop_know_fob = knows_other_dewormed / obs_know_person,
    prop_know_sob = thinks_other_knows / obs_know_person
  ) %>%
  left_join(
    original_distance_df,
    by = "original_cluster_id",
    relationship = "many-to-one"
  ) %>%
  mutate(
    assigned_treatment = assigned.treatment,
    assigned_dist_group = factor(original_dist_group, levels = c("close", "far")),
    signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
    signal = factor(signal, levels = c("no signal", "signal"))
  )

assert_true(!anyDuplicated(endline_know_A_itt_data$KEY.individ), "KEY.individ is duplicated in original-distance observability ITT data.")
assert_true(!any(is.na(endline_know_A_itt_data$assigned_dist_group)), "Missing original distance assignments for observability ITT clusters.")
assert_true(!any(is.na(endline_know_A_itt_data$mu_d)), "Missing expected-distance controls in original-distance observability ITT data.")
assert_true(
  endline_know_A_itt_data %>% distinct(original_cluster_id, cluster_id_rank) %>% count(original_cluster_id) %>% filter(n > 1) %>% nrow() == 0,
  "An observability original cluster ID maps to multiple bootstrap cluster ranks."
)

observability_source_check <- endline_know_A_itt_data %>%
  select(
    KEY.individ,
    original_cluster_id,
    assigned.treatment,
    dist.pot.group,
    dist.to.pot,
    standard_cluster.dist.to.pot
  ) %>%
  left_join(
    endline_data %>%
      transmute(
        KEY.individ,
        endline_cluster_id_raw = as.character(cluster.id),
        source_cluster_id = clean_cluster_id(cluster.id, "endline source cluster.id"),
        source_treatment = as.character(assigned.treatment),
        source_dist_group = as.character(dist.pot.group),
        source_dist_to_pot = dist.to.pot,
        source_standard_cluster_dist_to_pot = standard_cluster.dist.to.pot
      ),
    by = "KEY.individ",
    relationship = "one-to-one"
  )

assert_true(!any(is.na(observability_source_check$source_cluster_id)), "Some observability rows did not match endline source data.")
assert_true(all(observability_source_check$original_cluster_id == observability_source_check$source_cluster_id), "Observability cluster IDs differ from endline source data.")
assert_true(all(as.character(observability_source_check$assigned.treatment) == observability_source_check$source_treatment), "Observability treatment assignments differ from endline source data.")
assert_true(all(as.character(observability_source_check$dist.pot.group) == observability_source_check$source_dist_group), "Observability current distance groups differ from endline source data.")
assert_true(max(abs(observability_source_check$dist.to.pot - observability_source_check$source_dist_to_pot), na.rm = TRUE) < 1e-8, "Observability individual PoT distances differ from endline source data.")
assert_true(max(abs(observability_source_check$standard_cluster.dist.to.pot - observability_source_check$source_standard_cluster_dist_to_pot), na.rm = TRUE) < 1e-8, "Observability standardized cluster PoT distances differ from endline source data.")

observability_merge_checks <- tibble(
  check = c(
    "observability_rows",
    "observability_clusters",
    "missing_original_distance_assignments",
    "missing_expected_distance_controls",
    "treatment_mismatches_endline",
    "distance_group_mismatches_endline",
    "max_abs_dist_to_pot_diff_endline",
    "max_abs_standard_cluster_dist_to_pot_diff_endline"
  ),
  value = c(
    nrow(endline_know_A_itt_data),
    n_distinct(endline_know_A_itt_data$original_cluster_id),
    sum(is.na(endline_know_A_itt_data$assigned_dist_group)),
    sum(is.na(endline_know_A_itt_data$mu_d)),
    sum(as.character(observability_source_check$assigned.treatment) != observability_source_check$source_treatment),
    sum(as.character(observability_source_check$dist.pot.group) != observability_source_check$source_dist_group),
    max(abs(observability_source_check$dist.to.pot - observability_source_check$source_dist_to_pot), na.rm = TRUE),
    max(abs(observability_source_check$standard_cluster.dist.to.pot - observability_source_check$source_standard_cluster_dist_to_pot), na.rm = TRUE)
  )
)

write_csv(observability_merge_checks, "temp-data/original-distance-observability-merge-checks.csv")

original_distance_discrete_observability_regression <- function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] +  mu_d | county,
    data = data,
    weights = weights
  )
}

original_distance_fob_output <- wrapper_function(
  data = endline_know_A_itt_data %>%
    mutate(prop_knows = prop_know_fob),
  regression_spec = original_distance_discrete_observability_regression,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_original_distance_itt_fob_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/original-distance-itt-fob-tidy-tes.csv"
)

original_distance_sob_output <- wrapper_function(
  data = endline_know_A_itt_data %>%
    mutate(prop_knows = prop_know_sob),
  regression_spec = original_distance_discrete_observability_regression,
  table_options = list(
    dependent_var = "Dependent variable: Observability Beliefs"
  ),
  table_name = "rf_original_distance_itt_sob_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/original-distance-itt-sob-tidy-tes.csv"
)

original_distance_fob_output$tidy_summary %>%
  filter(assigned_dist_group == "far - close") %>%
  print(n = Inf)

original_distance_sob_output$tidy_summary %>%
  filter(assigned_dist_group == "far - close") %>%
  print(n = Inf)

}) # end original-distance ITT



#### NP Fit Plots -------------------------------------------------------------------

if (run_all || script_options$plots) run_section("NP Fit Plots", {

library(ggthemes)
canva_palette_vibrant = "Primary colors with a vibrant twist"
colour_tibble = tibble(
  assigned_treatment = c("control", "ink", "calendar", "bracelet", "signal"),
  colors = c(canva_pal(canva_palette_vibrant)(4), "#000000")
) %>%
  mutate(
    assigned_treatment = str_to_title(assigned_treatment),
    assigned_treatment = factor(assigned_treatment, levels = c("Control", "Ink", "Calendar", "Bracelet", "Signal"))
  )

summ_analysis_data = analysis_data %>%
  group_by(cluster.id, assigned.treatment) %>%
  summarise(
    mean_dist = mean(cluster.dist.to.pot),
    mean_dewormed = mean(dewormed)
  ) %>%
  ungroup() %>%
  left_join(colour_tibble, by = c("assigned.treatment" = "assigned_treatment")) %>%
  mutate(
    assigned.treatment = str_to_title(assigned.treatment),
    assigned.treatment = factor(assigned.treatment, levels = c("Control", "Ink", "Calendar", "Bracelet", "Signal"))
  )
max_dist = analysis_data %>%
  summarise(
    max_dist = max(cluster.dist.to.pot/1000)
  ) %>%
  pull(max_dist)

  analysis_data %>%
    group_by(assigned.treatment) %>%
    summarise(n_cl = n_distinct(cluster.id)) 
  analysis_data %>%
    group_by(assigned.treatment, dist.pot.group) %>%
    summarise(n_cl = n_distinct(cluster.id)) 

p_np_fit = cov_analysis_data %>% 
  mutate(
    assigned.treatment = str_to_title(assigned.treatment),
    assigned.treatment = factor(assigned.treatment, levels = c("Control", "Ink", "Calendar", "Bracelet", "Signal"))
  ) %>%
  ggplot(aes(
    x = dist.to.pot/1000,
    y = as.numeric(dewormed),
    colour = assigned.treatment
  )) +
  geom_smooth(se = FALSE, method = "gam") +
  geom_point(
    data = summ_analysis_data,
    aes(
      x = mean_dist/1000,
      y = mean_dewormed,
      colour = assigned.treatment
    )
  ) +
  theme_bw() +
  scale_colour_manual("", values = deframe(colour_tibble))  +
  scale_x_continuous(
    breaks = seq(0, 2.5, 0.5),
    labels = scales::number_format(accuracy = 0.1),
    limits = c(0, max(2.5, max_dist)) 
  ) +
  labs(
    x = "Distance to Treatment (km)",
    y = "Take-up Probability"
  ) +
  theme(legend.position = "bottom")
p_np_fit
ggsave(
  plot = p_np_fit,
  filename = file.path(
    "presentations",
    "takeup-np-rf-fit.pdf"
  ),
  width = 8,
  height = 6
)

}) # end plots



#### Change in Externality Knowledge (Table 13)

if (run_all || script_options$externality) run_section("Externality Knowledge", {

externality_data = endline_data %>%
    mutate(
      fully_aware_externalities = case_when(
        neighbours_worms_affect == "yes" & worms_affect == "yes" ~ TRUE, 
        # Ed: 2025-08-08 NA in these two variables is actually "don't know" due to 
        # a coding error in `analysis_util.R:129` in SurveyCTO these two 
        # variables use different binary encoding for yes/no and the original 
        # code corrects this but doesn't correct "don't know" correctly
        is.na(neighbours_worms_affect) | is.na(worms_affect) ~ FALSE,
        TRUE ~ FALSE
      ),
      know_worms_infectious = spread_worms == "yes",
      externality_omnibus = fully_aware_externalities | know_worms_infectious,
      assigned_treatment = assigned.treatment,
      assigned_dist_group = dist.pot.group
    ) 
    # select(KEY.individ, cluster.id, externality_omnibus) 


externality_knowledge_df = externality_data

full_externality_knowledge_df = full_analysis_data %>%
  mutate(
    female = gender == "female",
    cluster_id = dense_rank(cluster.id)
    ) %>%
  select(
    cluster_id,
    cluster.id,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group,
    # mu_d,
    # standard_cluster.dist.to.pot,
    # standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot),
    county,
    all_of(l_cov_vars),
    KEY.individ,
    everything()
    )  %>%
    inner_join(
      externality_data %>% select(-cluster.id),
      by = "KEY.individ"
    ) %>%
    left_join(
      cov_analysis_data %>%
        select(cluster.id.x, mu_d, standard_cluster.dist.to.pot) %>%
        unique(),
      by = c("cluster.id" = "cluster.id.x")
    )



externality_knowledge_regression = function(data, weights) {
  feols(
    externality_omnibus ~ 0 + assigned_treatment*assigned_dist_group + .[l_cov_vars] + mu_d   | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}


externality_knowledge_output = wrapper_function(
  data = externality_knowledge_df,
  regression_spec = externality_knowledge_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/externality-knowledge-tidy-tes.csv",
  table_name = "rf_externality_knowledge_tbl",
  table_options = list(
    dependent_var = "Dependent variable: Externality Knowledge"
  )
)


externality_knowledge_output$tidy_summary %>%
  select(
    assigned_treatment, assigned_dist_group, estimate, std_error, pval
  ) %>%
  print(n = 100)

#### Praise/Stigma by Externality Knowledge

praise_stigma_perception_data = baseline_data %>%
  select(cluster.id, matches("^(praise|stigma)_[^_]+$")) %>%
  gather(key = key, value = response, -cluster.id) %>%
  separate(key, c("praise.stigma", "topic"), "_") %>%
  separate(topic, c("topic", "question.group"), -2) %>%
  filter(!is.na(response))

cluster_perception_data = praise_stigma_perception_data %>%
  group_by(cluster.id, praise.stigma, topic) %>%
  summarise(
    mean_yes = mean(response == "yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    full_externality_knowledge_df %>%
      group_by(cluster.id) %>%
      summarise(
        mean_externality_knowledge = mean(externality_omnibus, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "cluster.id"
  )

praise_stigma_by_ext_know_fit = cluster_perception_data %>%
  filter(topic == "dewor") %>%
  ungroup() %>%
  mutate(above_median_externality_know = mean_externality_knowledge > median(mean_externality_knowledge, na.rm = TRUE)) %>%
  feols(
    mean_yes ~ above_median_externality_know,
    data = .,
    split = ~praise.stigma,
    vcov = "HC1"
  )

fitstat_register(
  "ctrl_mean_ps",
  function(est) {
    y = model.matrix(est, type = "lhs")
    list(mean = mean(y, na.rm = TRUE), sd = sprintf("(%.3f)", sd(y, na.rm = TRUE)))
  },
  alias = c("ctrl_mean_ps.mean" = "Dep. var. mean", "ctrl_mean_ps.sd" = ""),
  subtypes = c("mean", "sd")
)

etable(
  praise_stigma_by_ext_know_fit,
  tex = TRUE,
  title = "Praise and Stigma by Externality Knowledge",
  dict = c(
    mean_yes = "Fraction Yes",
    above_median_externality_knowTRUE = "Above Median Ext. Knowledge"
  ),
  fitstat = ~ctrl_mean_ps.mean + ctrl_mean_ps.sd + n,
  depvar = FALSE,
  digits = 3,
  digits.stats = 3,
  drop.section = "fixef",
  postprocess.tex = tex_postprocessing,
  replace = TRUE,
  style.df = style.df(depvar.title = "", fixef.title = "", var.title = "", stats.title = ""),
  notes = ""
  # file = file.path(
  #   params$table_output_path, "new-praise-stigma-by-externality-knowledge.tex"
  # )
)

}) # end externality


#### Takeup Continuous Distance + LASSO Covs + Cluster Expected Distance


# TODO: run balance on these people

if (run_all || script_options$takeup) run_section("Takeup Regressions", {

# First, attrition by treatment
cto_dropped_df = 
census_data %>%
  filter(!is.na(cluster.id)) %>%
  filter(sms.treatment == "sms.control") %>%
  filter(monitored) %>%
  filter(have_phone == "No" | have_phone == "Don't know number") %>%
  mutate(
    survey_cto_dropped = monitored & !true.monitored
  ) %>%
  filter(survey_cto_dropped)  %>%
  select(
    KEY.individ, cluster.id
  ) %>%
  left_join(
    analysis_data %>%
      transmute(
        cluster.id = as.numeric(levels(cluster.id)[cluster.id]), 
        county, assigned_treatment, assigned_dist_group,
        cluster_id, standard_cluster.dist.to.pot,
        cluster_id_rank
      ) %>%
      distinct(),
    by = "cluster.id"
  )

full_cto_dropped_df = bind_rows(
  cto_dropped_df %>% 
    mutate(dropped_cto = TRUE),
  analysis_data %>% mutate(cluster.id = as.numeric(levels(cluster.id))[cluster.id]) %>%
    mutate(dropped_cto = FALSE)
)



att_reg = function(data, weights) {
  feols(
    dropped_cto ~ 0 + assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

attrition_cto_output = wrapper_function(
  data = full_cto_dropped_df,
  regression_spec = att_reg,
  tidy_summ_path = "temp-data/tidy-rf-tes/attrition-cto-tidy-tes.csv",
  table_name = "rf_attrition_cto_tbl",
  table_options = list(
    dependent_var = "Dependent variable: Dropped from monitored sample"
  )
)



dist_cts_regression = function(data, weights) {
  feols(
    dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, cluster.dist.to.pot, "control")  +  mu_d + .[l_cov_vars] | county, 
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

dist_cts_output = wrapper_function(
  data = cov_analysis_data,
  regression_spec = dist_cts_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-dist-cts-tidy-tes.csv",
  table_name = "rf_dist_cts_spec_tbl"
)

dist_cts_exclude_baseline_output = wrapper_function(
  data = cov_analysis_data %>%
    filter(sms.treatment == "sms.control"),
  regression_spec = dist_cts_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-dist-cts-exclude-baseline-tidy-tes.csv",
  table_name = "rf_dist_cts_exclude_baseline_spec_tbl"
)

dist_cts_output$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)

dist_cts_exclude_baseline_output$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)

#### Takeup Continuous Distance + No Covs + Cluster Expected Distance
dist_cts_no_covs_regression = function(data, weights) {
  feols(
    dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, cluster.dist.to.pot, "control")  +  mu_d  | county, 
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

dist_cts_no_covs_output = wrapper_function(
  data = cov_analysis_data,
  regression_spec = dist_cts_no_covs_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-dist-cts-no-covs-tidy-tes.csv",
  table_name = "rf_dist_cts_no_covs_spec_tbl"
)



#### Takeup Discrete Distance + LASSO Covs + Cluster Expected Distance
discrete_distance_regression = function(data, weights) {
  feols(
    dewormed ~ 0 + assigned_treatment*assigned_dist_group + .[l_cov_vars] + mu_d | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

discrete_distance_covs_output = wrapper_function(
  data = cov_analysis_data,
  regression_spec = discrete_distance_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/discrete-dist-covs-tidy-tes.csv",
  table_name = "rf_discrete_dist_covs_tbl"
)

#### Takeup Discrete Distance + No LASSO Covs + No Cluster Expected Distance
discrete_distance_no_covs_no_mu_d_regression = function(data, weights) {
  feols(
    dewormed ~ 0 + assigned_treatment*assigned_dist_group | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

discrete_distance_no_covs_no_mu_d_output = wrapper_function(
  data = cov_analysis_data,
  regression_spec = discrete_distance_no_covs_no_mu_d_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/discrete-dist-no-covs-no-mu-d-tidy-tes.csv",
  table_name = "rf_discrete_dist_no_covs_no_mu_d_tbl"
)

##### Excluding baseline sample

discrete_distance_covs_exclude_baseline_output = wrapper_function(
  data = cov_analysis_data %>%
    filter(sms.treatment == "sms.control"),
  regression_spec = discrete_distance_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/discrete-dist-covs-exclude-baseline-tidy-tes.csv",
  table_name = "rf_discrete_dist_covs_exclude_baseline_tbl"
)

discrete_distance_covs_output$tidy_summary %>%
  print(n = Inf)

discrete_distance_covs_exclude_baseline_output$tidy_summary %>%
  print(n = Inf)


#### Takeup HH Distance + LASSO Covs + Cluster Expected Distance
hh_spec_regression = function(data, weights) {
  feols(
    dewormed ~  0  + assigned_treatment + dist.to.pot + i(assigned_treatment, dist.to.pot, "control")  + mu_d + .[l_cov_vars]  | county, 
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

hh_spec_output = wrapper_function(
  data = cov_analysis_data,
  regression_spec = hh_spec_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/hh-dist-tidy-tes.csv",
  table_name = "rf_hh_spec_tbl"
)

}) # end takeup regressions


#### Beliefs -------------------------------------------------------------------




if (run_all || script_options$beliefs) run_section("Beliefs Regressions", {

summ_know_A_df = summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A"))

summ_know_B_df = summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.B"))

endline_know_table_data %>%
  filter(fct_match(know.table.type, "table.A"))  %>%
  filter(num.recognized > 0) %>%
  group_by(
    relationship,
    dewormed
  ) %>%
  count() %>%
  pivot_wider(
    names_from = dewormed,
    values_from = n,
    values_fill = 0
  ) %>%
  rename(dont_know = `don't know`) %>%
  mutate(
    total = no + yes + dont_know
  ) %>%
  mutate(across(
    c(no, yes, dont_know),
    ~round(100*.x/total, 1)
  )) %>%
  arrange(dont_know) %>%
  mutate(
    relationship = case_when(
      relationship == "hh member" ~ "Household member",
      relationship == "extended family" ~ "Extended family",
      relationship == "friend" ~ "Friend",
      relationship == "neighbor" ~ "Neighbor",
      relationship == "church" ~ "Church member",
      relationship == "village member" ~ "Village member",
      relationship == "other" ~ "Other"
    )
  ) %>%
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    linesep = "",
    col.names = c(
      "Relationship",
      "No (\\%)",
      "Yes (\\%)",
      "Don't Know (\\%)",
      "Total"
    ),
    align = "lcccc"
  ) %>% 
  kable_styling(
    latex_options = c("scale_down")
  ) %>%
    custom_save_latex_table(
      table_name = "beliefs_relationships_tbl"
    )


belief_ana_df = analysis_data %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  filter(obs_know_person > 0)


belief_shared_vars = c(
  "KEY.individ", "assigned.treatment", "assigned_dist_group",
  "obs_know_person", "knows_other_dewormed_yes", "knows_other_dewormed_no",
  "doesnt_know_other_dewormed", "thinks_other_knows_yes", "thinks_other_knows_no",
  "doesnt_think_other_knows", "cluster.id", "cluster.dist.to.pot",
  "standard_cluster.dist.to.pot", "dist.to.pot", "county", "dewormed"
)


  
# Used for Conditional accuracy and Correct observability
endline_belief_df = endline_data %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  inner_join(summ_know_A_df, by = "KEY.individ") %>%
  select(any_of(belief_shared_vars), sms.treatment, dist.pot.group, assigned_treatment, assigned_dist_group,
         all_of(l_cov_vars), 
         standard_cluster.dist.to.pot,
         standard_clust_expected_dist, contains("know"), contains("pct"), mu_d, contains("cluster"))  %>%
  mutate(cluster_id = cluster.id)




# IDs in summ_know_A_df that aren't in endline_data
anti_join(summ_know_A_df, endline_data, by = "KEY.individ") %>%
  pull(KEY.individ)  %>% unique() %>% length()

# IDs in endline_data that aren't in summ_know_A_df
anti_join(endline_data, bind_rows(summ_know_A_df, summ_know_B_df), by = "KEY.individ") %>%
  pull(KEY.individ)  %>% unique() %>% length()


summ_know_A_df %>%
  filter(obs_know_person > 0) %>%
  summarize(
    n_indiv = n_distinct(KEY.individ)
  )

endline_know_A_df = endline_data %>%
  left_join(
    summ_know_A_df,
    by = "KEY.individ"
  ) %>%
  filter(sms.treatment == "sms.control", obs_know_person > 0) %>%
  mutate(
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group,
    prop_know_fob = knows_other_dewormed / obs_know_person,
    prop_know_sob = thinks_other_knows / obs_know_person
  ) 



discrete_pct_yesno = function(data, weights) {
  feols(
    pct_correct_classification_yesno ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] +  mu_d | county,
    data = data,
    weights = weights
  )
}

discrete_pct_yesnodk = function(data, weights) {
  feols(
    pct_correct_classification_yesnodk ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] +  mu_d | county,
    data = data,
    weights = weights
  )
}


discrete_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] +  mu_d | county,
    data = data,
    weights = weights
  )
}

discrete_f_know_no_covs_no_mu_d = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") | county,
    data = data,
    weights = weights
  )
}

cts_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") + .[l_cov_vars] + mu_d | county,
    data = data,
    weights = weights
  )
}

hh_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + dist.to.pot + i(assigned_treatment, dist.to.pot, "control")  + .[l_cov_vars] + mu_d | county,
    data = data,
    weights = weights
  )
}


#### FOB Discrete Distance + LASSO Covs + Cluster Expected Distance

# source("scratch/old-knowledge-reduced-form-code.R")


discrete_fob_output = wrapper_function(
  data = endline_know_A_df %>%
    mutate(prop_knows = prop_know_fob),
  regression_spec = discrete_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_discrete_fob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-fob-tidy-tes.csv"
)

#### FOB Discrete Distance + No LASSO Covs + No Cluster Expected Distance
discrete_fob_no_covs_no_mu_d_output = wrapper_function(
  data = endline_know_A_df %>%
    mutate(prop_knows = prop_know_fob),
  regression_spec = discrete_f_know_no_covs_no_mu_d,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_discrete_fob_no_covs_no_mu_d_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-fob-no-covs-no-mu-d-tidy-tes.csv"
)



# ── Lee bounds for discrete_fob_output ───────────────────────────────────────
# lee_base: full population for computing weighted selection rates per bootstrap draw.
# Column names aligned to assigned_treatment / assigned_dist_group to match endline_know_A_df.

lee_base = endline_data %>%
  filter(sms.treatment == "sms.control") %>%
  filter(!(KEY.individ %in% summ_know_B_df$KEY.individ)) %>%
  mutate(
    in_fob_sample       = KEY.individ %in% endline_know_A_df$KEY.individ,
    assigned_treatment  = assigned.treatment,
    assigned_dist_group = dist.pot.group
  )

fob_lee_upper = wrapper_function(
  data            = endline_know_A_df %>% mutate(prop_knows = prop_know_fob),
  regression_spec = discrete_f_know,
  lee_direction   = "upper",
  lee_base_data   = lee_base,
  tidy_summ_path  = "temp-data/tidy-rf-tes/fob-lee-upper-tidy-tes.csv",
  table_name      = "fob_lee_upper_tbl",
  table_options   = list(dependent_var = "Dependent variable: FOB (Lee upper bound)")
)

fob_lee_lower = wrapper_function(
  data            = endline_know_A_df %>% mutate(prop_knows = prop_know_fob),
  regression_spec = discrete_f_know,
  lee_direction   = "lower",
  lee_base_data   = lee_base,
  tidy_summ_path  = "temp-data/tidy-rf-tes/fob-lee-lower-tidy-tes.csv",
  table_name      = "fob_lee_lower_tbl",
  table_options   = list(dependent_var = "Dependent variable: FOB (Lee lower bound)")
)


fob_lee_lower$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)
fob_lee_upper$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)

discrete_fob_output$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)

# Build Lee bounds table: same layout as nice_kbl_table but showing [lb, ub] intervals
lee_bounds_tbl = function(lower_output, upper_output, caption, dependent_var) {
  tbl_dist_levels     = c("combined", "close", "far", "far - close")
  tbl_contrast_levels = c("bracelet", "calendar", "ink", "control", "Observations")

  lower_est = lower_output$tidy_summary %>%
    select(assigned_treatment, assigned_dist_group, lb = estimate, se_lb = std_error, pval, n_obs_line, show_pval_only)

  upper_est = upper_output$tidy_summary %>%
    select(assigned_treatment, assigned_dist_group, ub = estimate, se_ub = std_error)

  # For close/far/combined: lb from lower model, ub from upper model
  cell_bounds = lower_est %>%
    left_join(upper_est, by = c("assigned_treatment", "assigned_dist_group")) %>%
    filter(assigned_dist_group != "far - close")

  # Far-close: correct interval arithmetic
  # lb(far - close) = lb(far) - ub(close)
  # ub(far - close) = ub(far) - lb(close)
  # SEs combined in quadrature
  far_close_bounds = lower_est %>%
    filter(assigned_dist_group %in% c("far", "close")) %>%
    left_join(upper_est, by = c("assigned_treatment", "assigned_dist_group")) %>%
    select(assigned_treatment, assigned_dist_group, lb, ub, se_lb, se_ub, pval, n_obs_line, show_pval_only) %>%
    pivot_wider(names_from = assigned_dist_group, values_from = c(lb, ub, se_lb, se_ub, pval)) %>%
    mutate(
      lb      = lb_far - ub_close,
      ub      = ub_far - lb_close,
      se_lb   = sqrt(se_lb_far^2 + se_ub_close^2),
      se_ub   = sqrt(se_ub_far^2 + se_lb_close^2),
      pval    = pval_far,
      assigned_dist_group = "far - close"
    ) %>%
    select(assigned_treatment, assigned_dist_group, lb, ub, se_lb, se_ub, pval, n_obs_line, show_pval_only)

  tbl = bind_rows(cell_bounds, far_close_bounds) %>%
    mutate(
      estim_std = case_when(
        n_obs_line ~ pval,
        TRUE       ~ linebreak(
          paste0("[", round(lb, 3), ", ", round(ub, 3), "]",
                 "\n", "(", round(se_lb, 3), "), (", round(se_ub, 3), ")"),
          align = "c")
      ),
      assigned_dist_group = factor(assigned_dist_group, tbl_dist_levels),
      assigned_dist_group = fct_relabel(assigned_dist_group, str_to_title),
      assigned_treatment  = factor(assigned_treatment, tbl_contrast_levels)
    ) %>%
    filter(assigned_treatment %in% tbl_contrast_levels) %>%
    arrange(assigned_dist_group, assigned_treatment) %>%
    select(assigned_treatment, assigned_dist_group, estim_std) %>%
    pivot_wider(names_from = assigned_dist_group, values_from = estim_std) %>%
    mutate(assigned_treatment = fct_relabel(assigned_treatment, str_to_title))

  tbl %>%
    kbl(
      col.names = c(dependent_var, paste0("(", 1:4, ")")),
      format    = "latex",
      linesep   = "\\addlinespace",
      booktabs  = TRUE,
      escape    = FALSE,
      align     = "lcccc",
      caption   = caption
    ) %>%
    kable_styling(latex_options = "scale_down") %>%
    add_header_above(c(" ", "Combined", "Close", "Far", "Far - Close"), line = FALSE) %>%
    add_header_above(c(" " = 1, "Lee Bounds" = 4)) %>%
    row_spec(3, hline_after = TRUE)
}

fob_lee_tbl = lee_bounds_tbl(
  lower_output  = fob_lee_lower,
  upper_output  = fob_lee_upper,
  caption       = "Lee Bounds: First-Order Beliefs",
  dependent_var = "Dependent variable: Observability"
)

fob_lee_tbl %>%
  custom_save_latex_table(table_name = "fob_lee_bounds_tbl")
stop()
# ─────────────────────────────────────────────────────────────────────────────

discrete_pct_yesno_output = wrapper_function(
  data = endline_belief_df,
  regression_spec = discrete_pct_yesno,
  table_options = list(
    dependent_var = "Dependent variable: Conditional accuracy (Yes/No)"
  ),
  table_name = "rf_discrete_pct_yesno_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-pct-yesno-tidy-tes.csv"
)


discrete_pct_yesnodk_output = wrapper_function(
  data = endline_belief_df,
  regression_spec = discrete_pct_yesnodk,
  table_options = list(
    dependent_var = "Dependent variable: Conditional accuracy (Yes/No/Don't Know)"
  ),
  table_name = "rf_discrete_pct_yesnodk_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-pct-yesnodk-tidy-tes.csv"
)


panel_a_tex_path = file.path(
  script_options$table_output_path,
  "rf_discrete_fob_spec_tbl.tex"
)

panel_b_tex_path = file.path(
  script_options$table_output_path,
  "rf_discrete_pct_yesno_spec_tbl.tex"
)

panel_c_tex_path = file.path(
  script_options$table_output_path,
  "rf_discrete_pct_yesnodk_spec_tbl.tex"
)

#' Combine multiple .tex panel tables into a single table with Panel A/B/C headers.
#'
#' @param panel_paths Named character vector of .tex file paths. Names become
#'   panel titles (e.g., "Panel A: Observability").
#' @param output_path Path to write the combined .tex file.
#' @param notes Optional character string of table notes (placed inside a
#'   threeparttable tablenotes environment). NULL for no notes.
#' @return Invisibly returns the combined tex lines.
combine_panel_tables <- function(panel_paths, output_path, notes = NULL,
                                 drop_rows = NULL, replacements = NULL) {

  # Read all panels
  panels <- lapply(panel_paths, readLines)

  # Helper: extract body lines between first \midrule and \bottomrule
  extract_body <- function(lines) {
    midrule_idx <- which(grepl("^\\\\midrule", lines))
    bottomrule_idx <- which(grepl("^\\\\bottomrule", lines))
    stopifnot(length(midrule_idx) >= 1, length(bottomrule_idx) >= 1)
    # Body is everything from the first \midrule to just before \bottomrule
    lines[midrule_idx[1]:max(bottomrule_idx - 1)]
  }

  # Count columns from the first panel's tabular spec
  first_panel <- panels[[1]]
  tabular_line <- first_panel[grep("\\\\begin\\{tabular\\}", first_panel)[1]]
  # Count 'c' and 'l' in the column spec to get number of columns
  col_spec <- regmatches(tabular_line, regexpr("\\{[lcr]+\\}", tabular_line))
  n_cols <- nchar(gsub("[^lcr]", "", col_spec))

  # Build the header from the first panel (everything up to and excluding the
  # first \midrule), but remove the "Dependent variable:" row
  header_end <- which(grepl("^\\\\midrule", first_panel))[1] - 1
  header_lines <- first_panel[1:header_end]
  header_lines <- header_lines[!grepl("^Dependent variable:", header_lines)]

  # Build combined output
  out <- header_lines

  for (i in seq_along(panels)) {
    panel_title <- names(panel_paths)[i]
    body <- extract_body(panels[[i]])
    # Drop unwanted rows (e.g., |Calendar| - |Bracelet|)
    if (!is.null(drop_rows)) {
      body <- body[!grepl(drop_rows, body)]
    }
    # Skip only the first \midrule (we add our own before the panel title),
    # but keep internal \midrules (e.g., between coefficients and statistics)
    first_midrule <- which(grepl("^\\\\midrule$", body))[1]
    if (!is.na(first_midrule)) {
      body <- body[-first_midrule]
    }
    # Apply text replacements
    if (!is.null(replacements)) {
      for (j in seq_along(replacements)) {
        body <- gsub(names(replacements)[j], replacements[[j]], body)
      }
    }
    # Panel title row spanning all columns
    out <- c(out,
      sprintf("\\midrule"),
      sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", n_cols, panel_title),
      body
    )
  }

  out <- c(out, "\\bottomrule")

  # Close tabular
  out <- c(out, "\\end{tabular}}")

  # Wrap in threeparttable if notes are provided
  if (!is.null(notes)) {
    out <- c(
      "\\begin{threeparttable}",
      out,
      "\\begin{tablenotes}[flushleft]",
      "\\small",
      sprintf("\\item \\textit{Notes:} %s", notes),
      "\\end{tablenotes}",
      "\\end{threeparttable}"
    )
  }

  writeLines(out, output_path)
  message("Combined panel table written to: ", output_path)
  invisible(out)
}

# Combine the three belief decomposition panels
combine_panel_tables(
  panel_paths = c(
    "Panel A: Definite-answer rate (Table 1 observability)" = panel_a_tex_path,
    "Panel B: Conditional accuracy (given definite answer)" = panel_b_tex_path,
    "Panel C: Correct observability" = panel_c_tex_path
  ),
  output_path = file.path(
    script_options$table_output_path,
    "rf_discrete_belief_decomposition_tbl.tex"
  ),
  drop_rows = "\\$\\|Calendar\\|",
  replacements = c("^Control " = "Control mean ")
)



#### FOB Continuous Distance + LASSO Covs + Cluster Expected Distance
cts_fob_output = wrapper_function(
  data = endline_know_A_df %>%
    mutate(prop_knows = prop_know_fob),
  regression_spec = cts_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_cts_fob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-cts-fob-tidy-tes.csv"
)

#### FOB HH Distance + LASSO Covs + Cluster Expected Distance
hh_fob_output = wrapper_function(
  data = endline_know_A_df %>%
    mutate(prop_knows = prop_know_fob),
  regression_spec = hh_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_hh_fob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-hh-fob-tidy-tes.csv"
)

#### SOB Discrete Distance + LASSO Covs + Cluster Expected Distance
discrete_sob_output = wrapper_function(
  data = endline_know_A_df %>%
    mutate(prop_knows = prop_know_sob),
  regression_spec = discrete_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability Beliefs"
  ),
  table_name = "rf_discrete_sob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-sob-tidy-tes.csv"
)

}) # end beliefs regressions



#### Levels --------------------------------------------------------------------

#### Takeup LEVELS Continuous Distance + LASSO Covs + Expected Distance

if (run_all || script_options$takeup) run_section("Takeup Levels", {

# For plotting
dist_cts_spec_bs_draws = map_dfr(
  1:500,
  ~bayes_bs_f(
    seed = .x,
    f = dist_cts_regression,
    data = cov_analysis_data
  ),
  .progress = TRUE
  )

dist_cts_spec_levels = actual_bayesian_bs_fit(
  seed = "realised fit",
  f = dist_cts_regression,
  data = cov_analysis_data
)$bs_fit %>%
  filter(!is.na(assigned_treatment)) 

dist_cts_spec_levels_ci = dist_cts_spec_bs_draws %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  summarise(
    conf.low = quantile(mean_pred, 0.025),
    conf.high = quantile(mean_pred, 0.975)
  ) %>%
  filter(!is.na(assigned_treatment))

tidy_dist_cts_spec_levels = left_join(
  dist_cts_spec_levels,
  dist_cts_spec_levels_ci,
  by = c("assigned_treatment", "assigned_dist_group")
) %>%
  select(-signal, -seed) %>%
  rename(estimate = mean_pred)
tidy_dist_cts_spec_levels %>%
  write_csv("temp-data/reducedform-dist-cts-tidy-levels.csv")  


#### Takeup LEVELS Discrete Distance + LASSO Covs + Expected Distance
discrete_distance_covs_bs_draws = map_dfr(
  1:500,
  ~bayes_bs_f(
    seed = .x,
    f = discrete_distance_regression,
    data = cov_analysis_data
  ),
  .progress = TRUE
  )

discrete_distance_covs_levels = actual_bayesian_bs_fit(
  seed = "realised fit",
  f = discrete_distance_regression,
  data = cov_analysis_data
)$bs_fit %>%
  filter(!is.na(assigned_treatment)) 

discrete_distance_covs_levels_ci = discrete_distance_covs_bs_draws %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  summarise(
    conf.low = quantile(mean_pred, 0.025),
    conf.high = quantile(mean_pred, 0.975)
  ) %>%
  filter(!is.na(assigned_treatment))

tidy_discrete_distance_cov_levels = left_join(
  discrete_distance_covs_levels,
  discrete_distance_covs_levels_ci,
  by = c("assigned_treatment", "assigned_dist_group")
) %>%
  select(-signal, -seed) %>%
  rename(estimate = mean_pred)
tidy_discrete_distance_cov_levels %>%
  write_csv("temp-data/discrete-dist-covs-tidy-levels.csv")

}) # end takeup levels



#### FOB Levels Discrete Distance + LASSO Covs + Expected Distance

if (run_all || script_options$beliefs) run_section("Beliefs Levels", {

fob_bs_draws = map_dfr(
  1:500,
  ~bayes_bs_f(
    seed = .x,
    f = discrete_f_know,
    data = endline_know_A_df %>%
      mutate(prop_knows = prop_know_fob)
  ),
  .progress = TRUE
  )

fob_levels_point = actual_bayesian_bs_fit(
  seed = "realised fit",
  f = discrete_f_know,
  data = endline_know_A_df %>%
    mutate(prop_knows = prop_know_fob)
)$bs_fit %>%
  filter(!is.na(assigned_treatment)) 

fob_levels_ci = fob_bs_draws %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  summarise(
    conf.low = quantile(mean_pred, 0.025),
    conf.high = quantile(mean_pred, 0.975)
  ) %>%
  filter(!is.na(assigned_treatment))

fob_levels = left_join(
  fob_levels_point,
  fob_levels_ci,
  by = c("assigned_treatment", "assigned_dist_group")
) %>%
  select(-signal, -seed) %>%
  rename(estimate = mean_pred)

fob_levels %>%
  write_csv("temp-data/reducedformfob-tidy-levels.csv")

}) # end beliefs levels



#### Alternative Regressions ---------------------------------------------------

if (run_all || script_options$endline) run_section("Endline/Incentive/Preference/Travel", {

####  Endline Predicted Deworming Takeup
endline_data = endline_data %>%
  mutate(
    assigned_treatment = as_factor(assigned.treatment), 
    assigned_dist_group = as_factor(dist.pot.group),
    cluster_id = as_factor(cluster.id),
    # this isn't actually used
    standard_cluster.dist.to.pot = dist.to.pot/sd_of_dist,
    dworm_frac = dworm_rate / 10,
    # different naming convention here
    have_ink = ink_visible
  ) 

mis_endline_df = endline_data %>%
left_join(
  summ_endline_know_table %>% transmute(KEY.individ, know.table.type, obs_know_person),
  by = "KEY.individ"
) %>%
  # filter(obs_know_person > 0 | is.na(obs_know_person))  %>% 
  filter(know.table.type == "table.A" | is.na(know.table.type)) %>%
  mutate(attrit_know_table = is.na(know.table.type))


mis_endline_df %>%
  count(attrit_know_table)

summ_know_A_df = summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A"))

endline_data %>%
  inner_join(
summ_know_A_df %>%
  filter(obs_know_person > 0) , by = "KEY.individ"

  )



# Knowledge table attrition
endline_data = endline_data %>%
  mutate(
    attrit_know_table = !in_know_table
  )
know_table_fit = function(data, weights) {
  feols(
    attrit_know_table ~ 0 + assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

know_table_A_fit = function(data, weights) {
  feols(
    attrit_know_table ~ 0 + assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + is_backup | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

mis_endline_df = mis_endline_df %>%
  mutate(is_backup = KEY.individ %in% census_data$KEY.individ[census_data$endline.backup == TRUE])

know_table_attrit_output = wrapper_function(
  data = endline_data,
  regression_spec = know_table_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/know-table-attrition-tidy-tes.csv",
  table_name = "know_table_attrition_spec_tbl",
  table_options = list(
    caption = "Knowledge Table Attrition", 
    dependent_var = "Dependent variable: Not in knowledge table"
    )
)

know_table_attrit_output$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)

know_table_A_attrit_output = wrapper_function(
  data = mis_endline_df,
  regression_spec = know_table_A_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/know-table-A-attrition-tidy-tes.csv",
  table_name = "know_table_A_attrition_spec_tbl",
  table_options = list(
    caption = "Knowledge Table Attrition",
    dependent_var = "Dependent variable: Not in knowledge table"
    )
)

know_table_A_attrit_output$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)


pred_dworm_fit = function(data, weights) {
  feols(
    dworm_frac ~ 0 + assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] + mu_d | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}


pred_dworm_output = wrapper_function(
  data = endline_data,
  regression_spec = pred_dworm_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/predicted-endline-deworm-takeup-tidy-tes.csv",
  table_name = "predicted_endline_deworm_takeup_spec_tbl",
  table_options = list(
    caption = "Average Treatment Effects: Reduced Form", 
    dependent_var = "Dependent variable: Predicted Take-up", 
    type = "APE", 
    stars = TRUE,
    drop_H0s = TRUE
    )
)

#### Incentive Implementation --------------------------------------------------

mean_deworm_string_f = function(string) {
  str_detect(str_to_lower(string), "drug|medicine|tablet|deworm|Deworm|worm|treat|dewom|deform|medic")
}

endline_data = endline_data %>%
  mutate(
    meandeworm_bracelet = mean_deworm_string_f(bracelet_meaning),
    meandeworm_ink = mean_deworm_string_f(ink_meaning),
    meandeworm_cal = mean_deworm_string_f(cal_meaning)
  )



bad_vars = c(
  "bad_bracelet",
  "bad_ink",
  "bad_calendar"
)


got_vars = c(
  "got_bracelet", 
  "got_ink", 
  "got_cal"
)
have_vars = c(
  "have_bracelet", 
  "have_cal", 
  "have_ink"
)
seen_vars = c(
  "seen_bracelet", 
  "seen_ink", 
  "seen_cal"
)
mean_vars = c(
  "meandeworm_bracelet", 
  "meandeworm_ink", 
  "meandeworm_cal"
)

display_vars = c(
  "display_bracelet",
  "display_ink",
  "display_cal"
)

endline_data %>%
  mutate(
    have_incentive = coalesce(have_bracelet, have_ink, have_cal),
    have_incentive_phone = coalesce(have_bracelet, have_ink, have_cal, !is.na(have_phone))
  ) %>%
  count(have_incentive_phone)

endline_data %>%
  mutate(
    have_incentive = coalesce(have_bracelet, have_ink, have_cal),
    have_incentive_phone = coalesce(have_bracelet, have_ink, have_cal, have_phone == "Yes")
  ) %>%
  count(have_incentive)

endline_data = endline_data %>%
  mutate(
    display_bracelet = wear_bracelet == 1,
    display_ink = ink_visible == 1,
    display_cal = case_when(
      str_detect(
        str_to_lower(cal_where), 
        "wall|living|table|hang"
      ) ~ TRUE,
      !is.na(cal_where) ~ FALSE,
      TRUE ~ NA
    )
  )


long_incentive_check_df = endline_data %>%
  select(all_of(c(got_vars, have_vars, seen_vars, mean_vars, display_vars, bad_vars)), assigned_treatment, cluster_id, county)  %>%
  pivot_longer(
    cols = all_of(c(got_vars, have_vars, seen_vars, mean_vars, display_vars, bad_vars))
  ) %>%
  mutate(
    variable_type = str_extract(name, "(\\w+)(?=_)"),
    name = str_extract(name, "(?<=_)\\w+"), 
    name = if_else(name == "cal", "calendar", name)
    )   %>%
  filter(name == assigned_treatment)  %>%
  mutate(
    treat_type = paste0(assigned_treatment, "_", variable_type)
  ) %>%
  select(-name)


long_incentive_check_df %>%
  group_by(variable_type) %>%
  summarise(
    n_na = sum(is.na(value)),
    n_non_na = sum(!is.na(value))
  )

tidy_incentive_check_df = long_incentive_check_df %>%
  filter(!is.na(value)) %>%
  feols(
    value ~ i(assigned_treatment, "ink") ,
    split = ~variable_type,
    cluster = ~cluster_id
  ) %>%
  map_dfr(
    ~tidy(.x) %>%
    mutate(n = nobs(.x)), 
    .id = "lhs"
  ) %>%
  mutate(
    treatment = str_extract(
      term, "(?<=assigned_treatment::)\\w+"
    ),
    treatment = replace_na(treatment, "ink")
  ) %>%
  mutate(
    variable_type = str_extract(
      lhs, 
      "(?<=sample: )\\w+$"
    )
    ) %>%
  select(
    -lhs,
    -term
  )



wide_incentive_check_input_df = tidy_incentive_check_df %>%
  mutate(across(c(estimate, std.error), round, 3)) %>%
  mutate(
    stars = case_when(
      treatment != "ink" & p.value < 0.01 ~ "***",
      treatment != "ink" & p.value < 0.05 ~ "**", 
      treatment != "ink" & p.value < 0.1 ~ "*" , 
      TRUE ~ ""
    ),
    estim_std = linebreak(paste0(estimate, stars, "\n", str_glue("({std.error})")), align = "c") 
  ) %>%
  mutate(treatment = str_replace(treatment, "ink", "ink (levels)")) %>%
  select( 
    treatment, 
    variable_type, 
    estim_std
  ) %>%
  mutate(treatment = str_to_title(treatment)) %>%
  spread(
    variable_type, 
    estim_std
  ) %>%
  select(
    treatment,
    seen, # "have you seen people wearing these "X"?"/ "have you seen people these calendars before"
    meandeworm, # "What does it mean if a person has a bracelet?" -> coded into deworm mentions
    got, # "Did you receive X when you went for deworming?"
    have, # "Do you still have X?"
    display, # "Do you wear the bracelet?" / "Is the ink visible?"
    bad # Drawbacks of incentive
  )




wide_incentive_check_input_df = wide_incentive_check_input_df %>%
  bind_rows(
    long_incentive_check_df %>%
      filter(!is.na(value))  %>%
      group_by(
        variable_type
      ) %>%
      summarise(
        n = n()
      ) %>%
      mutate(n = as.character(prettyNum(n, big.mark = ","))) %>%
      pivot_wider(
        names_from = variable_type,
        values_from = n
      ) %>%
      mutate(
        treatment = "Observations"
      )
  )

incentive_check_tbl = wide_incentive_check_input_df %>%
  knitr::kable(
    format = "latex",
      col.names = c(
        "", 
        "Seen incentive",
        "Link to deworming",
        "Received incentive",
        "Have incentive",
        "Wearing/displayed",
        "Reported any drawback"
        ),
      escape = FALSE, 
      booktabs = TRUE,
      align = "lccccc", 
      caption = "Endline Incentive Checks"
  ) %>% 
  row_spec(c(2), hline_after = TRUE) 


incentive_check_tbl %>%
  custom_save_latex_table("incentive-check-tbl")

bracelet_distance_data = endline_data %>%
  filter(assigned_treatment == "bracelet") %>%
  mutate(
    assigned_dist_group = factor(assigned_dist_group, levels = c("close", "far")),
    seen_in_community = seen_bracelet == 1,
    received_bracelet = dewormed.reported == 1 & got_bracelet == 1,
    still_has_bracelet = if_else(
      received_bracelet,
      have_bracelet == 1,
      NA
    ),
    wearing_bracelet = if_else(
      received_bracelet,
      wear_bracelet == 1,
      NA
    )
  ) %>%
  select(
    cluster_id,
    assigned_dist_group,
    seen_in_community,
    still_has_bracelet,
    wearing_bracelet
  ) %>%
  pivot_longer(
    cols = c(seen_in_community, still_has_bracelet, wearing_bracelet),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c("seen_in_community", "still_has_bracelet", "wearing_bracelet"),
      labels = c("Seen in community", "Still has bracelet", "Wearing bracelet")
    )
  )

bracelet_distance_level_fit = feols(
  value ~ 0 + assigned_dist_group,
  split = ~outcome,
  data = bracelet_distance_data,
  cluster = ~cluster_id
)

bracelet_distance_diff_fit = feols(
  value ~ i(assigned_dist_group, ref = "close"),
  split = ~outcome,
  data = bracelet_distance_data,
  cluster = ~cluster_id
)

bracelet_distance_levels = bracelet_distance_level_fit %>%
  map_dfr(
    ~tidy(.x) %>%
      select(term, estimate, std.error, p.value),
    .id = "outcome"
  ) %>%
  mutate(
    dist_group = case_when(
      term == "assigned_dist_groupclose" ~ "Close",
      term == "assigned_dist_groupfar" ~ "Far"
    )
  )

bracelet_distance_diff = bracelet_distance_diff_fit %>%
  map_dfr(
    ~tidy(.x) %>%
      filter(term == "assigned_dist_group::far") %>%
      select(estimate, std.error, p.value),
    .id = "outcome"
  ) %>%
  mutate(dist_group = "Far--Close")

bracelet_distance_tbl_input = bind_rows(
  bracelet_distance_levels,
  bracelet_distance_diff
) %>%
  mutate(
    stars = case_when(
      dist_group != "Far--Close" ~ "",
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      p.value < 0.1 ~ "*",
      TRUE ~ ""
    ),
    estimate_se = linebreak(
      paste0(round(estimate, 3), stars, "\n", "(", round(std.error, 3), ")"),
      align = "c"
    ),
    dist_group = factor(dist_group, levels = c("Close", "Far", "Far--Close"))
  ) %>%
  select(outcome, dist_group, estimate_se) %>%
  pivot_wider(
    names_from = dist_group,
    values_from = estimate_se
  ) %>%
  mutate(
    outcome = str_remove(outcome, "^sample\\.var: outcome; sample: "),
    outcome = factor(
      outcome,
      levels = c(
        "Seen in community",
        "Still has bracelet",
        "Wearing bracelet"
      ),
      labels = c("Seen in community", "Still has bracelet", "Wearing bracelet")
    )
  ) %>%
  arrange(outcome)

bracelet_distance_tbl = bracelet_distance_tbl_input %>%
  knitr::kable(
    format = "latex",
    col.names = c("", "Close", "Far", "Far--Close"),
    escape = FALSE,
    booktabs = TRUE,
    align = "lccc",
    caption = "Bracelet Visibility, Retention, and Wearing by Distance\\label{tab:bracelet-distance-display}"
  ) %>%
  kable_styling(
    latex_options = c("scale_down")
  )

bracelet_distance_tbl %>%
  custom_save_latex_table("bracelet-distance-display-tbl")

bracelet_distance_tbl_path = file.path(
  params$table_output_path,
  "bracelet-distance-display-tbl.tex"
)
bracelet_distance_tbl_tex = readLines(bracelet_distance_tbl_path, warn = FALSE)
bracelet_distance_tbl_tex = c(
  "\\caption{Bracelet Visibility, Retention, and Wearing by Distance}",
  "\\label{tab:bracelet-distance-display}",
  bracelet_distance_tbl_tex,
  "\\floatfoot{\\textit{Notes:} The table reports endline measures of bracelet visibility, possession, and wearing among respondents assigned to the Bracelet arm. Outcomes are reported separately for Close and Far communities, along with Far--Close differences. ``Seen in community'' equals one if the respondent reports having seen the bracelet in the community. ``Still has bracelet'' and ``Wearing bracelet'' are defined among respondents who dewormed and report receiving the bracelet. Standard errors are clustered at the community level.}"
)
writeLines(bracelet_distance_tbl_tex, bracelet_distance_tbl_path)

# incentive drawbacks
incentive_drawback_tbl = wide_incentive_check_input_df %>%
  select(treatment, bad) %>%
  knitr::kable(
    format = "latex",
    col.names = c(
      "Treatment", 
      "Reported any drawback"
    ),
    escape = FALSE, 
    booktabs = TRUE,
    align = "lc", 
    caption = "Endline Incentive Drawbacks"
  ) %>% 
  row_spec(c(2), hline_after = TRUE) 

incentive_drawback_tbl %>%
  custom_save_latex_table("incentive-drawback-tbl")
 

long_incentive_check_df %>%
  group_by(variable_type) %>%
  summarise(
    n_na = sum(is.na(value)),
    n_non_na = sum(!is.na(value))
  )
### Preference for Gift Fit Not Dewormed ---------------------------------------

pref_gift_fit_not_dewormed_full_sample = full_analysis_data %>%
    # 38,019
    filter(!is.na(gift_choice)) %>%
    # 3,676
    filter(monitored) %>%
    # 3,329
    filter(monitor.consent) %>%
    # 3,329 %>%
    filter(dewormed == FALSE) %>%
    # 1,808
    filter(sms.treatment == "sms.control") %>%
    # 1,365
    filter(gift_choice %in% c("bracelet", "calendar")) %>%
    # 1,350
    group_by(assigned.treatment, dist.pot.group, dewormed) %>% 
    mutate(arm.size = n()) %>% 
    group_by(gift_choice, .add = TRUE) %>%
    ungroup() %>%
    select(KEY.individ, cluster.id, gift_choice, 
      assigned.treatment, dist.pot.group, county, gender,
      age.census
      ) %>%
    mutate(
      want_bracelet = gift_choice == "bracelet"
    )  %>%
    mutate(
      assigned_treatment = factor(assigned.treatment),
      assigned_dist_group = factor(dist.pot.group),
      cluster_id = factor(cluster.id),
      female = gender == "female"
      ) %>%
      left_join(
        cov_analysis_data %>%
          select(cluster.id.x, mu_d, standard_cluster.dist.to.pot) %>%unique(),
          by = c("cluster.id" = "cluster.id.x")
      ) %>%
      ungroup() %>%
      mutate(cluster_id_rank = dense_rank(cluster.id))


summ_gift_df = full_analysis_data %>%
    filter(!is.na(gift_choice)) %>%
    # 3,676
    filter(monitored) %>%
    # 3,329
    filter(monitor.consent) %>%
    # 3,329 %>%
    filter(dewormed == FALSE) %>%
    # 1,808
    group_by(assigned.treatment, dist.pot.group, dewormed) %>% 
    mutate(arm.size = n()) %>% 
    group_by(gift_choice, add = TRUE) %>%
    ungroup() %>%
    select(KEY.individ, cluster.id, gift_choice, 
      assigned.treatment, dist.pot.group, county, gender,
      age.census
      ) %>%
    mutate(
      want_bracelet = gift_choice == "bracelet"
    )  %>%
    mutate(assigned_treatment = factor(assigned.treatment))
    

summ_gift_df  %>%
    group_by(assigned_treatment, dist.pot.group) %>%
    summarise(
      n = n(),
      mean_want_bracelet = mean(want_bracelet)
    ) %>%
    bind_rows(
      summ_gift_df %>%
        group_by(assigned_treatment) %>%
        summarise(
          dist.pot.group = "combined",
          n = n(),
          mean_want_bracelet = mean(want_bracelet)
        )
    ) %>%
    arrange(assigned_treatment) %>%
    group_by(dist.pot.group) %>%
    mutate(
      te = mean_want_bracelet - mean_want_bracelet[assigned_treatment == "control"], 
    )



pref_gift_fit_not_dewormed = analysis_data %>%
    filter(
      !is.na(gift_choice), 
      monitored, 
      monitor.consent, 
      !hh.baseline.sample.pool, 
      !is.na(sms.treatment)) %>% 
    group_by(assigned.treatment, dist.pot.group, dewormed) %>% 
    mutate(arm.size = n()) %>% 
    group_by(gift_choice, add = TRUE) %>%
    filter(
      dewormed == FALSE
    )  %>%
    ungroup() %>%
    select(
      KEY.individ, cluster.id, gift_choice, assigned.treatment, dist.pot.group,
      county, standard_cluster.dist.to.pot, cluster_id_rank) %>%
    mutate(
      want_bracelet = gift_choice == "bracelet"
    )  %>%
    mutate(
      assigned_treatment = factor(assigned.treatment),
      assigned_dist_group = factor(dist.pot.group),
      cluster_id = factor(cluster.id)
      ) %>%
      left_join(
        cov_analysis_data %>%
          select(KEY.individ, mu_d, all_of(l_cov_vars)),
          by = "KEY.individ"
      )



pref_gift_fit = function(data, weights) {
  feols(
    want_bracelet ~  assigned_treatment*assigned_dist_group  + .[l_cov_vars] + mu_d | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
} 


pref_fit_flipped = wrapper_function(
  data = pref_gift_fit_not_dewormed,
  regression_spec = pref_gift_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/preference-for-bracelet-FLIPPED-tidy-tes.csv",
  table_name = "preference_for_bracelet_FLIPPED_spec_tbl",
  table_options = list(
    caption = "Average Treatment Effects: Reduced Form", 
    dependent_var = "Dependent variable: Prefer Bracelet", 
    stars = TRUE
    ),
    flip_calendar_sign = TRUE
)


pref_fit_not_flipped = wrapper_function(
  data = pref_gift_fit_not_dewormed,
  regression_spec = pref_gift_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/preference-for-bracelet-tidy-tes.csv",
  table_name = "preference_for_bracelet_spec_tbl",
  table_options = list(
    caption = "Average Treatment Effects: Reduced Form", 
    dependent_var = "Dependent variable: Prefer Bracelet", 
    stars = TRUE
    )
)
  pref_fit_not_flipped$tidy_summary %>%       
    filter(                               
      assigned_treatment %in% c("bracelet", "calendar", "ink", "control", "abs(calendar) - abs(bracelet)") |
      n_obs_line                                                                                               
    ) %>%                                 
    prep_tbl(stat = params$stat, stars = TRUE) %>%                                                             
    mutate(                                                                                                    
      assigned_treatment = str_replace(                                                                      
        as.character(assigned_treatment),                                                                      
        fixed("$|Calendar| - |Bracelet|$"),                                                                  
        "$|\\text{Calendar}| - |\\text{Bracelet}|$"                                                          
      )
    ) %>%                                                                                                      
    arrange(factor(assigned_treatment, levels = c(
      "Bracelet", "Calendar", "Ink", "Control",                                                                
      "$|\\text{Calendar}| - |\\text{Bracelet}|$",                                                           
      "Observations"                                                                                         
    ))) %>%                                                                                                    
    nice_kbl_table(
      cap = "...",                                                                                             
      outcome_var = "Dependent variable: Prefer Bracelet"                                                    
    ) %>%                                                                                                    
    custom_save_latex_table(table_name = "preference_for_bracelet_abs_diff_spec_tbl") 
pref_fit_not_flipped$tidy_summary %>%
  filter(str_detect(assigned_treatment, "cal|bra")) %>%
  print(n = 100)


pref_not_flipped_full = wrapper_function(
  data = pref_gift_fit_not_dewormed_full_sample,
  regression_spec = pref_gift_fit,
  tidy_summ_path = "temp-data/temp.csv",
  table_name = "temp-pref-bra",
  table_options = list(
    caption = "Average Treatment Effects: Reduced Form", 
    dependent_var = "Dependent variable: Prefer Bracelet", 
    stars = TRUE
    ),
    flip_calendar_sign = FALSE
)

pref_not_flipped_full

pref_not_flipped_full$tidy_summary %>%
  filter(str_detect(assigned_treatment, "cal|bra")) %>%
  print(n = 100)



pref_gift_not_dewormed_full_sample_fit = wrapper_function(
  data = pref_gift_fit_not_dewormed_full_sample,
  regression_spec = pref_gift_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/preference-for-bracelet-full-sample-FLIPPED-tidy-tes.csv",
  table_name = "preference_for_bracelet_full_sample_FLIPPED_spec_tbl",
  table_options = list(
    caption = "Average Treatment Effects: Reduced Form", 
    dependent_var = "Dependent variable: Prefer Bracelet", 
    stars = TRUE
    ),
    flip_calendar_sign = TRUE
)

pref_gift_not_dewormed_full_sample_fit$tidy_summary %>%
  print(n = 100)



#### Travel Time --------------------------------------------------------

endline_data %>%
  count(travel)

  endline_data %>%
  mutate(
    travel_clean = case_when(
      travel == "1" ~ "walk",
      travel == "2" ~ "motorbike",
      travel == "3" ~ "car/taxi",
      travel == "4" ~ "bus",
      travel == "5" ~ "free ride"
    )
  ) %>%
    group_by(dist.pot.group) %>%
    summarise(
      pr_walk = mean(travel == "1", na.rm = TRUE),
      ed = mean(travel_clean == "walk", na.rm = TRUE),
      more_robust_pr_walk = mean(str_detect(travel, "1"), na.rm = TRUE)
    )


travel_time_df = endline_data %>%
  mutate(
    travel = if_else(travel == "99", NA_character_, travel),
  ) %>%
  bind_rows(
    .,
      mutate(., dist.pot.group = "combined")
  ) %>%
  group_by(dist.pot.group) %>%
  summarise(
    pr_walk = mean(str_detect(travel, "1"), na.rm = TRUE),
    mean_time_travel = mean(time_travel, na.rm = TRUE),
    median_time_travel = median(time_travel, na.rm = TRUE),
    sd_time_travel = sd(time_travel, na.rm = TRUE),
    q_lo = quantile(time_travel, 0.25, na.rm = TRUE),
    q_hi = quantile(time_travel, 0.75, na.rm = TRUE),
    n = sum(!is.na(time_travel))
  ) %>%
  pivot_longer(
    -dist.pot.group,
    names_to = "stat",
    values_to = "value"
  )  %>%
  mutate(
    value = case_when(
      stat %in% c("mean_time_travel", "median_time_travel", "sd_time_travel", "q_lo", "q_hi") ~ as.character(round(value, 0)),
      TRUE ~ as.character(prettyNum(round(value, 3), big.mark = ","))
    )
  ) %>%
  pivot_wider(
    names_from = dist.pot.group,
    values_from = value
  )   %>%
  bind_rows(
    analysis_data %>%
      bind_rows(
        .,
        mutate(., dist.pot.group = "combined")
      ) %>%
      group_by(dist.pot.group) %>%
      summarise(
        dist_mean = prettyNum(round(mean(dist.to.pot/1000), digits = 3), big.mark = ","),
      ) %>%
      pivot_wider(
        names_from = dist.pot.group,
        values_from = dist_mean
      ) %>%
      mutate(
        stat = "Distance to treatment (km)"
      ),
      .
  ) %>%
  mutate(
    stat  = case_when(
      stat == "mean_time_travel" ~ "Mean",
      stat == "median_time_travel" ~ "Median",
      stat == "sd_time_travel" ~ "Standard deviation",
      stat == "q_lo" ~ "25th Percentile",
      stat == "q_hi" ~ "75th Percentile",
      stat == "n" ~ "Observations",
      stat == "pr_walk" ~ "Fraction walked to treatment",
      stat == "Distance to treatment (km)" ~ "Distance to treatment (km)"
    )
  ) %>%
  select(stat, combined, close, far)


travel_time_df  %>%
  knitr::kable(
    format = "latex",
    col.names = c("", "Combined", "Close", "Far"),
    escape = FALSE, 
    booktabs = TRUE,
    align = "lcccccc", 
    caption = "Travel to Treatment Location",
    digits = 3
  ) %>%
  pack_rows(
    "Travel time (minutes)",
    3,8,
    italic = TRUE,
    escape = FALSE,
    # latex_gap_space = latex_group_gap_space, 
    hline_after = TRUE, 
    hline_before = TRUE,
    bold = TRUE
  ) %>%
  row_spec(c(7), hline_after = TRUE)  %>%
  custom_save_latex_table("travel-time-tbl")

analysis_data %>%
  unique() %>%
  group_by(assigned_dist_group) %>%
  summarise(
    n_in_close = mean(dist.to.pot < 1250, na.rm = TRUE),
    n_in_far = mean(dist.to.pot >= 1250, na.rm = TRUE)
  )

endline_data = endline_data %>%
  mutate(
    travel_clean = case_when(
      travel == "1" ~ "walk",
      travel == "2" ~ "motorbike",
      travel == "3" ~ "car/taxi",
      travel == "4" ~ "bus",
      travel == "5" ~ "free ride"
    )
  )


endline_data %>%
  count(travel_clean) %>%
  mutate(
    pct = 100*n/sum(n)
  )

endline_data %>%
  summarise(
    mean_time_travel = mean(time_travel, na.rm = TRUE),
    median_time_travel = median(time_travel, na.rm = TRUE)
  )


endline_data %>%
  summarise(
    frac_pay_0 = mean(travel_pay == 0, na.rm = TRUE),
    mean_pay = mean(travel_pay, na.rm = TRUE)
  )

endline_data %>%
  mutate(
    travel = case_when(
      travel == "1" ~ "walk",
      travel == "2" ~ "motorbike",
      travel == "3" ~ "car/taxi",
      travel == "4" ~ "bus",
      travel == "5" ~ "free ride"
    )
  ) %>%
  group_by(dist.pot.group) %>%
  count(travel) %>%
  arrange(-n) %>%
  pivot_wider(
    names_from = dist.pot.group,
    values_from = n
  ) %>%
  mutate(
    pct_close = 100*close / sum(close, na.rm = TRUE),
    pct_far = 100*far / sum(far, na.rm = TRUE)
  )

}) # end endline






#### SMS -----------------------------------------------------------------------

if (run_all || script_options$sms) run_section("SMS Heterogeneity", {
  stop()

pval_only_terms = c("bracelet - calendar", "signal")

sms_analysis_data %>%
  count(gender)

anti_join(
  cluster_expected_dist_df %>%
    select(cluster.id) %>%
    unique(),
  sms_analysis_data %>%
    select(cluster.id) %>%
    unique(),
  by = "cluster.id"
)

anti_join(
  sms_analysis_data %>%
    select(cluster.id) %>%
    unique(),
  cluster_expected_dist_df %>%
    select(cluster.id) %>%
    unique(),
  by = "cluster.id"
)

sms_analysis_data_control = sms_analysis_data %>%
  filter(assigned_treatment == "control") %>%
  left_join(
    cluster_expected_dist_df %>%
      transmute(cluster.id, mu_d = as.numeric(clust_expected_dist)/sd_of_dist),
    by = c("cluster.id" = "cluster.id")
  ) %>%
  mutate(
    sms_treatment = factor(sms_treatment, levels = c("smscontrol", "reminderonly", "socialinfo")),
    female = gender == "female",
    cluster.id = as.numeric(cluster.id)
    ) 
  
sms_control_reg = sms_analysis_data_control %>%
  feols(
    dewormed ~  sms_treatment + .[l_cov_vars] + mu_d | county,
    data = .,
    cluster = ~cluster_id
  )

fitstat_register(
  "control_mean_sms_control",
  function(est) {
    ctrl = sms_analysis_data_control %>% filter(sms_treatment == "smscontrol")
    ctrl_mean_se = feols(
      dewormed ~ 1,
      data = ctrl,
      cluster = ~cluster_id
    ) %>%
      coeftable() %>%
      .[1, "Std. Error"]
    list(mean = mean(ctrl$dewormed, na.rm = TRUE), se = sprintf("(%.3f)", ctrl_mean_se))
  },
  alias = c("control_mean_sms_control.mean" = "Control mean", "control_mean_sms_control.se" = ""),
  subtypes = c("mean", "se")
)

etable(
  list("Takeup" = sms_control_reg),
  tex = TRUE,
  title = "Effects of Reminder and Social Info SMS on Take-up in Control Group",
  dict = c(
    dewormed = "Take-up",
    sms_treatmentsocialinfo = "Social Info",
    sms_treatmentreminderonly = "Reminder Only"
  ),
  keep = c("Social Info", "Reminder Only"),
  fitstat = ~control_mean_sms_control.mean + control_mean_sms_control.se + n,
  depvar = FALSE,
  digits = 3,
  digits.stats = 3,
  drop.section = "fixef",
  postprocess.tex = compact_fixest_postprocessing,
  replace = TRUE,
  style.df = style.df(depvar.title = "", fixef.title = ""),
  notes = "",
  file = file.path(
    params$table_output_path, "incentive-control-sms-te.tex"
  )
)
# \begingroup
# \centering
# \begin{tabular}{lc}
#    \tabularnewline \midrule \midrule
#    Model:        & (1)\\  
#    \midrule
#    \emph{Variables}\\
#    Reminder Only & 0.167$^{***}$\\   
#                  & (0.029)\\   
#    Social Info   & 0.131$^{***}$\\   
#                  & (0.028)\\   
#    \midrule
#    \emph{Fit statistics}\\
#    Control mean  & 0.330\\  
#                  & (0.032)\\
#    Observations  & 1,705\\  
#    \midrule \midrule
#    \multicolumn{2}{l}{\emph{Clustered (cluster\_id) standard-errors in parentheses}}\\
#    \multicolumn{2}{l}{\emph{Signif. Codes: ***: 0.01, **: 0.05, *: 0.1}}\\
# \end{tabular}
# \par\endgroup





    # prop_knows ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] +  mu_d | county,
f_sms = function(data, weights) {
  feols(
    dewormed ~ 0  + 
      assigned_treatment + 
      standard_cluster.dist.to.pot + 
      sms_treatment + 
      i(assigned_treatment, standard_cluster.dist.to.pot, "control") +
      i(assigned_treatment, sms_treatment, "control") +
      i(sms_treatment, standard_cluster.dist.to.pot) +
      sms_treatment:assigned_treatment:standard_cluster.dist.to.pot 
      | county,
    data = data,
    weights = weights
  )
}


sms_bs_draws = map_dfr(
    1:500,
    ~bayes_bs_f(
        seed = .x, 
        f = f_sms, 
        data = sms_analysis_data,
        sms_treatment
    ),
    .progress = TRUE
)


clean_bs_sms_signal_draws = sms_bs_draws %>%
  clean_signal_draws(sms_treatment)

clean_bs_sms_te_draws = sms_bs_draws %>%
  clean_te_draws(sms_treatment)


create_sms_te = function(draws) {
  draws %>%
    group_by(seed, assigned_dist_group, sms_treatment) %>%
    mutate(
      te = if_else(assigned_treatment == "control", mean_pred, mean_pred - mean_pred[assigned_treatment == "control"])
    )  %>%
    ungroup() %>%
    group_by(seed, assigned_dist_group, assigned_treatment) %>%
    mutate(
      diff_te = te - te[sms_treatment == "smscontrol"]
    ) 
}

sms_bs_tes = sms_bs_draws %>%
  filter(!is.na(assigned_treatment)) %>%
  select(-signal) %>%
  add_predictions(sms_treatment)  %>%
  create_sms_te() %>%
  rename(estimate = diff_te)

sms_signal_bs_tes = sms_bs_draws %>%
  filter(!is.na(signal)) %>%
  select(-assigned_treatment) %>%
  add_signal_predictions(sms_treatment) %>%
  group_by(seed, assigned_dist_group, sms_treatment) %>%
    mutate(
      te = if_else(signal == "no signal", mean_pred, mean_pred - mean_pred[signal == "no signal"])
    )  %>%
    ungroup() %>%
    group_by(seed, assigned_dist_group, signal) %>%
    mutate(
      diff_te = te - te[sms_treatment == "smscontrol"]
   )  %>%
  rename(estimate = diff_te)  %>%
  rename(assigned_treatment = signal)


realised_sms_fit = actual_bayesian_bs_fit(
  seed = "realised fit",
  f = f_sms,
  data = sms_analysis_data,
  sms_treatment
)$bs_fit


realised_sms_tes = realised_sms_fit %>%
  filter(!is.na(assigned_treatment)) %>%
  select(-signal) %>%
  add_predictions(sms_treatment)  %>%
  create_sms_te() %>%
  ungroup() %>%
  rename(realised_pred = diff_te) %>%
  select(assigned_dist_group, assigned_treatment, sms_treatment, realised_pred)

realised_sms_signal_fit = realised_sms_fit %>%
  filter(!is.na(signal)) %>%
  select(-assigned_treatment) %>%
  add_signal_predictions(sms_treatment) %>%
  group_by(seed, assigned_dist_group, sms_treatment) %>%
    mutate(
      te = if_else(signal == "no signal", mean_pred, mean_pred - mean_pred[signal == "no signal"])
    )  %>%
    ungroup() %>%
    group_by(seed, assigned_dist_group, signal) %>%
    mutate(
      diff_te = te - te[sms_treatment == "smscontrol"]
   )  %>%
  rename(realised_pred = diff_te) %>%
  ungroup() %>%
  select(assigned_dist_group, assigned_treatment = signal, sms_treatment, realised_pred)


realised_sms_tes
realised_sms_signal_fit



both_sms_fits = bind_rows(
  sms_bs_tes,
  sms_signal_bs_tes
) %>%
  mutate(
    show_pval_only = assigned_treatment %in% pval_only_terms
  ) %>%
  filter(assigned_treatment != "no signal") 

realised_sms_both = bind_rows(
  realised_sms_signal_fit,
  realised_sms_tes
) 



    clean_sms_tes = both_sms_fits %>%
      group_by(
          assigned_treatment,
          assigned_dist_group,
          sms_treatment
      ) %>%
      summarise(
          std_error = sd(estimate),
          conf.low = quantile(estimate, (1 - ci_width)/2),
          conf.high = quantile(estimate, 1 - (1 - ci_width)/2)
      ) %>%
      left_join(
          realised_sms_both,
          by = c("assigned_dist_group", "assigned_treatment", "sms_treatment")
      ) %>%
      mutate(
          pval = 2*pnorm(-abs(realised_pred)/std_error),
          oneside_pval = pnorm(-realised_pred/std_error)
      ) %>%
      mutate(
          pval = round(pval, 4),
          oneside_pval = round(oneside_pval, 4)
      ) %>%
      select(
          assigned_treatment, 
          assigned_dist_group, 
          sms_treatment,
          realised_pred, 
          std_error, 
          conf.low,
          conf.high,
          pval, 
          oneside_pval) %>%
      rename(estimate = realised_pred)  %>%
      filter(sms_treatment != "smscontrol")

clean_sms_tes %>%
  write_csv("temp-data/differential-tes-by-sms.csv")

clean_sms_tes %>%
  filter(assigned_treatment != "control") %>%
  select(assigned_treatment, assigned_dist_group, sms_treatment, pval, oneside_pval)


# clean_sms_tes %>%
#   filter(sms_treatment != "smscontrol")  %>%
#   mutate(show_pval_only = FALSE) %>%
#   filter(sms_treatment != "reminderonly") %>%
#   mutate(
#     show_pval_only = assigned_treatment %in% pval_only_terms
#   ) %>%
#   prep_tbl(stat = params$stat) %>%
#   nice_kbl_table(
#     cap = "Heterogeneous SMS Average Treatment Effects",
#     outcome_var = "Dependent variable: Take-up"
#   ) %>%
#   custom_save_latex_table(
#     table_name = "sms_diff_tes_tbl"
#   )

library(ggthemes)


p_sms_tes = clean_sms_tes %>%
  filter(sms_treatment != "smscontrol")  %>%
  mutate(show_pval_only = FALSE)  %>%
  filter(assigned_treatment != "signal") %>%
  select(
    assigned_treatment,
    assigned_dist_group,
    sms_treatment,
    estimate,
    conf.low,
    conf.high
  ) %>%
  mutate(
    assigned_treatment = case_when(
      assigned_treatment == "bracelet - calendar" ~ "Bracelet - Calendar",
      assigned_treatment == "bracelet" ~ "Bracelet",
      assigned_treatment == "calendar" ~ "Calendar",
      assigned_treatment == "ink" ~ "Ink",
      assigned_treatment == "control" ~ "Control Mean",
    ),
    assigned_treatment = factor(
      assigned_treatment,
      levels = c(
        "Control Mean",
        "Bracelet - Calendar",
        "Ink",
        "Calendar",
        "Bracelet"
      )
    ),
    assigned_dist_group = str_to_title(assigned_dist_group),
    sms_treatment = case_when(
      sms_treatment == "smscontrol" ~ "SMS Control",
      sms_treatment == "reminderonly" ~ "Reminder Only",
      sms_treatment == "socialinfo" ~ "Social Info"
    )
  ) %>%
  ggplot(aes(
    x = estimate,
    xmin = conf.low,
    xmax = conf.high,
    y = assigned_treatment,
    colour = sms_treatment
  )) +
  geom_pointrange(
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~assigned_dist_group) +
  geom_vline(
    xintercept = 0,
    linetype = "longdash"
  ) +
  labs(
    x = "Estimate",
    y = "",
    colour = ""
  ) +
  scale_colour_canva(
    "",
    palette = "Primary colors with a vibrant twist"
  )

p_sms_tes
ggsave("temp-data/p-sms-tes.pdf", width = 8, height = 6)

}) # end sms

#### Heterogeneity -------------------------------------------------------------

if (run_all || script_options$heterogeneity) run_section("Heterogeneity + WTP", {

library(marginaleffects)
analysis_data = analysis_data %>%
  mutate(cluster.id = as.character(cluster.id)) %>%
  left_join(
    baseline_worm %>%
      mutate(cluster.id = as.character(cluster.id)) %>%
      group_by(cluster.id) %>%
      summarise(
        frac_externality = mean(fully_aware_externalities, na.rm = TRUE),
        frac_prev_dewormed = mean(treated_lgl, na.rm = TRUE),
      ) %>%
      mutate(
        frac_externality_gt_mean = frac_externality > mean(frac_externality, na.rm = TRUE),
        frac_prev_dewormed_gt_mean = frac_prev_dewormed > mean(frac_prev_dewormed, na.rm = TRUE),
        cluster.id = factor(cluster.id)
      ),
      by = "cluster.id"
  ) %>%
  mutate(
    age_gt_40 = age.census > 40,
    treatment = factor(assigned_treatment, levels = c("control", "bracelet", "calendar", "ink"))
  ) 


clean_perception_data = baseline_data %>% 
  select(cluster.id, matches("^(praise|stigma)_[^_]+$")) %>% 
  gather(key = key, value = response, -cluster.id) %>% 
  separate(key, c("praise.stigma", "topic"), "_") %>% 
  separate(topic, c("topic", "question.group"), -2)  %>%
  filter(!is.na(response))  


overall_judgement_score_df = clean_perception_data %>%
  count(cluster.id, praise.stigma, topic, response)  %>%
  group_by(cluster.id, praise.stigma, topic) %>% 
  mutate(n = n/sum(n))   %>%
  group_by(cluster.id) %>%
  filter(response == "yes") %>%
  summarise(
    judge_score = prod(n), 
  ) %>% 
  ungroup() %>%
  mutate(
    mean_judge_score = mean(judge_score, na.rm = TRUE)
  ) %>%
  mutate(topic = "overall") %>%
  mutate(
    judge_score_gt_mean = judge_score > mean_judge_score
  ) %>%
  mutate(cluster.id = factor(cluster.id))


het_ols = function(data, judge_data) {
  age_fit = data %>%
    feols(
      dewormed ~  
        treatment + 
        standard_cluster.dist.to.pot + 
        age_gt_40 +
        i(treatment, age_gt_40, "control")   +
        factor(county),
        cluster = ~cluster.id
    ) 
  judge_fit = data %>%
    left_join(
      judge_data %>% 
        select(cluster.id, judge_score_gt_mean),
        by = "cluster.id"
    ) %>%
    feols(
      dewormed ~  
        treatment + 
        standard_cluster.dist.to.pot + 
        judge_score_gt_mean +
        i(treatment, judge_score_gt_mean, "control") +
        factor(county),
        cluster = ~cluster.id
    )
  phone_fit = data %>%
    feols(
      dewormed ~  
        treatment + 
        standard_cluster.dist.to.pot + 
        have_phone_lgl +
        i(treatment, have_phone_lgl, "control") +
        factor(county),
        cluster = ~cluster.id
    ) 

  prevdeworm_fit = data %>%
    feols(
      dewormed ~  
        treatment + 
        standard_cluster.dist.to.pot + 
        frac_prev_dewormed_gt_mean +
        i(treatment, frac_prev_dewormed_gt_mean, "control") +
        factor(county),
        cluster = ~cluster.id
    )
    
  externality_fit = data %>%
    feols(
      dewormed ~  
        treatment + 
        standard_cluster.dist.to.pot + 
        frac_externality_gt_mean +
        i(treatment, frac_externality_gt_mean, "control")   +
        factor(county),
        cluster = ~cluster.id
    )


  gender_fit = data %>%
    mutate(female = gender == "female") %>%
    feols(
      dewormed ~  
        treatment + 
        standard_cluster.dist.to.pot + 
        female +
        i(treatment, female, "control") +
        factor(county),
        cluster = ~cluster.id
    )
    return(list(
      age_fit = age_fit,
      judge_fit = judge_fit,
      phone_fit = phone_fit,
      prevdeworm_fit = prevdeworm_fit,
      externality_fit = externality_fit,
      gender_fit = gender_fit
    ))
}  



het_fits = het_ols(analysis_data, overall_judgement_score_df)

age_het_fit = het_fits$age_fit
judge_het_fit = het_fits$judge_fit
phone_het_fit = het_fits$phone_fit
prevdeworm_het_fit = het_fits$prevdeworm_fit
externality_het_fit = het_fits$externality_fit
gender_het_fit = het_fits$gender_fit


  etable(
    list(
      "Age $>40$" = age_het_fit,
      "Female" = gender_het_fit,
      "Phone Owner" = phone_het_fit,
      "Judgement" = judge_het_fit,
      "Previously Dewormed" = prevdeworm_het_fit,
      "Externality Knowledge" = externality_het_fit
    ),
    drop = "dist|Constant|factor\\(county\\)",
    dict = c(
      frac_externality_gt_mean = "Covariate",
      frac_externality_gt_meanTRUE = "Covariate",

      age_gt_40 = "Covariate",
      age_gt_40TRUE = "Covariate",

      judge_score_gt_mean = "Covariate",
      judge_score_gt_meanTRUE = "Covariate",


      frac_prev_dewormed_gt_mean = "Covariate",
      frac_prev_dewormed_gt_meanTRUE = "Covariate",

      have_phone_lgl = "Covariate",
      have_phone_lglTRUE = "Covariate",


      female = "Covariate",
      femaleTRUE = "Covariate",


      treatmentink = "Ink",
      treatmentcalendar = "Calendar",
      treatmentbracelet = "Bracelet",
      "treatment = ink" = "Ink",
      "treatment = calendar" = "Calendar",
      "treatment = bracelet" = "Bracelet",

      "bracelet" = "Bracelet",
      "calendar" = "Calendar",
      "ink" = "Ink",
      "cluster.id" = "Community"
    ),
    fitstat = c("my", "n"),
    headers = list(
      "Age $>40$" = 1,
      "Female" = 1,
      "Phone Owner" = 1,
      "Judgemental" = 1,
      "Prev Dewormed" = 1,
      "Externality Knowledge" = 1
    ),
    depvar = FALSE,
    digits = 3,
    digits.stats = 3,
    file = file.path(
      params$table_output_path, "het-tes.tex"
    ),
    drop.section = "fixef",
    tex = TRUE, 
    postprocess.tex = compact_fixest_postprocessing,
    replace = TRUE,
    style.df = style.df(
      depvar.title = "", 
      fixef.title = "")
  )

# WTP Checks

wtp_choice_data <- wtp.data %>%
  select(-any_of("county")) %>%
  left_join(
    cluster.strat.data %>% select(cluster.id, county),
    by = "cluster.id"
  )
if (!"second_choice" %in% names(wtp_choice_data)) {
  wtp_choice_data <- wtp_choice_data %>%
    mutate(
      first_choice = factor(first_choice, levels = 1:2, labels = c("bracelet", "calendar")),
      second_choice = if_else(first_choice == "bracelet", cal_plus_ksh, bra_plus_ksh) %>%
        factor(levels = 1:2, labels = c("switch", "keep"))
    )
}

wtp_out = full_analysis_data %>%
  mutate(stratum = county) %>%
  prepare_bayes_wtp_data(
    wtp_choice_data,
    
    preference_value_diff = seq(-100, 100, 10), 
    num_preference_value_diff = length(preference_value_diff), 
    
    wtp_utility_df = 3,
    tau_mu_wtp_diff = 100,
    mu_wtp_df_student_t = 7,
    tau_sigma_wtp_diff = 50,
    sigma_wtp_df_student_t = 2.5
  )

 full_analysis_data %>%
    filter(
           assigned.treatment == "control") %>%
    filter(!is.na(gift_choice)) %>%
    mutate(
      sms_control = sms.treatment.2 == "sms.control",
      monitored = true.monitored
    ) %>%
    group_by(sms_control, monitored, gift_choice) %>%
    summarise(
      n = n()
    ) %>%
    pivot_wider(
      names_from = gift_choice,
      values_from = n
    )  %>%
    mutate(
      n_total = bracelet + calendar,
      pct_calendar = 100 * calendar / (calendar + bracelet),
      n_neither = neither
    ) %>%
    select(-neither)

 full_analysis_data %>%
    count(
      wtp_samp = !is.na(gift_choice), 
      sms_control = sms.treatment.2 == "sms.control",
      monitored = true.monitored
      ) 

origin_prepared_analysis_data = full_analysis_data
origin_prepared_analysis_data %>%
    filter(!is.na(gift_choice) ,
     gift_choice != "neither",
           assigned.treatment == "control",
           sms.treatment.2 == "sms.control")

  incentive_choice_data <- origin_prepared_analysis_data %>%
    filter(!is.na(gift_choice) & gift_choice != "neither",
           assigned.treatment == "control",
           sms.treatment.2 == "sms.control") %>%
    select(county, cluster.id, gift_choice, phone_owner) %>%
    mutate(offer = 0,
           response = "keep")

  incentive_choice_data <- wtp_choice_data %>%
    semi_join(origin_prepared_analysis_data, "cluster.id") %>% # Make sure we have the same clusters
    filter(!is.na(first_choice)) %>%
    transmute(county, cluster.id,
              gift_choice = first_choice,
              offer = price,
              response = second_choice) %>%
    bind_rows(incentive_choice_data) %>%
    mutate(gift_choice = 2 * (gift_choice == "calendar") - 1,
           response = 2 * (response == "switch") - 1) 

incentive_choice_data


wtp_choice_data %>%
  select(first_choice, price)

}) # end heterogeneity
