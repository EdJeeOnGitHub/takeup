library(tidyverse)
library(posterior)
library(tidybayes)
library(broom)
library(data.table)
library(knitr)
library(kableExtra)
library(ggthemes)
library(fixest)
library(magrittr)
library(stringr)

if (interactive()) {
    params = lst(
        table_output_path = "presentations/rf-tables",
        show_probs = FALSE,
        width = 0.95,
        cache = FALSE,
        fit = FALSE,
        stat = "std.error" # "ci", "p", "std.error"
    )
    source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
    source(file.path("analysis_util.R"))
    source(file.path( "dist_structural_util.R"))
    source(file.path("multilvlr", "multilvlr_util.R"))
} else {
    source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
    source(file.path("analysis_util.R"))
    source(file.path( "dist_structural_util.R"))
    source(file.path("multilvlr", "multilvlr_util.R"))
}

models_we_want = c(params$struct_models, params$rf_models)

options(
    dplyr.show_progress = FALSE, 
    digits = 4, 
    knitr.kable.NA = '')

knitr::opts_chunk$set(
    echo = FALSE, 
    cache = params$cache, 
    warnings = FALSE,
    warning = FALSE,
    message = FALSE,
    cache.path = stringr::str_glue("rf-takeup-table-cache/"), 
    fig.path = str_glue("rf-takeup-table-fig-cache/"), 
    fig.align = "center")

dir.create("presentations/rf-tables", showWarnings = FALSE)


ci_width = as.numeric(params$width)
treat_levels_c = c("control", "ink", "calendar", "bracelet")
treat_levels = c("ink", "calendar", "bracelet")
dist_levels = c("close", "far")
model_level_order = c("reduced form", "structural")

quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

output_basepath = file.path(
  params$output_path,
  str_glue("output_dist_fit{params$fit_version}")
)


canva_palette_vibrant <- "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))
load(file.path("data", "takeup_village_pot_dist.RData"))
load(file.path("data", "analysis.RData"))

baseline.data = read_rds("temp-data/reclean_baseline_data.rds") # Not sampling data!

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()



monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


analysis_data <- monitored_nosms_data

## Load Census Data
census_data_env = new.env()
with_env = function(f, e = parent.frame()) {
    stopifnot(is.function(f))
    environment(f) = e
    f
}
load_census_function = function(){
  load(file.path("data", "takeup_census.RData"))
  return(census.data)
}
census_data = with_env(load_census_function, census_data_env)() %>%
  rename(census.consent = consent) # Rename this to reduce chance of error

n_indiv_df = census_data %>%
    group_by(cluster.id) %>%
    summarise(
        n_per_cluster = sum(num.individuals)
    )


clean_worm_covariates = function(data) {
  cov_data = data %>%
    unnest(when_treat) %>%
    mutate(
      who_worms_other = if_else(is.na(who_worms_other), "", who_worms_other),
      stop_worms_other = replace_na(stop_worms_other, ""),
      when_treat = replace_na(when_treat, "DK"),
      # these are nested lists of responses so we map_lgl and use any()
      all_can_get_worms = map2_lgl(
        who_worms, 
        who_worms_other,
        ~any(
          "everyone" %in% .x | 
          str_detect(str_to_lower(.y), "any") | 
          (("adult" %in% .x) & ("child" %in% .x)) |
          ("adult" %in% .x | str_detect(str_to_lower(.y), "adult|man|woman|men|women|person")) & ("child" %in% .x | str_detect(str_to_lower(.y), "child|under|young|teenager|below"))
      )
      ), 
      correct_when_treat = when_treat %in% c("every 3 months",
                                     "every 6 months"), 
      know_how_stop_worms = map2_lgl(
        stop_worms, 
        stop_worms_other,
        ~any(.x %in% c(
          "medicine", 
          "wearing shoes", 
          "using toilets", 
          "wash hands") | str_detect(.y, "cooked|prepar|cook"))),
      adult_in_family_treated = who_treated %in% c("adult", "both")
    )

  cov_data = cov_data %>%
    mutate(
      fully_aware_externalities = case_when(
        neighbours_worms_affect == "yes" & worms_affect == "yes" ~ TRUE, 
        is.na(neighbours_worms_affect) | is.na(worms_affect) ~ NA,
        TRUE ~ FALSE
      ),
      know_worms_infectious = spread_worms == "yes"
    )
  
    treated_past_present = "treated" %in% colnames(data) 
    if (treated_past_present) {
      cov_data = cov_data %>%
        mutate(
          treated_lgl = case_when(
            treated == "yes" ~ TRUE, 
            treated == "no" ~ FALSE, 
            TRUE ~ NA
          ),
      family_treated_lgl = family_treated == "yes"
        )
    }
    return(cov_data)
}

clean_pretreat_covariates = function(baseline_data, endline_data) {
  drop_vars = c(
    "starttime", 
    "endtime", 
    "age_census", 
    "phone_census", 
    "return",
    "no_return",
    "interview",
    "no_interview",
    "language_return",
    "isValidated"
    )
  cov_data = bind_rows(
    # baseline_data %>% select(-any_of(drop_vars)) %>% mutate(sample_type = "baseline"),
    endline_data %>% select(-any_of(drop_vars)) %>% mutate(sample_type = "endline")
  ) 
  cov_data = cov_data %>%
    mutate(
      floor_tile_cement = floor == "Cement" | floor == "Tiles"
    )
  school_year_df = tribble(
    ~school, ~years_schooling, 
    "Never gone to school", 0,
    "Primary 1", 1,
    "Primary 2", 2,
    "Primary 3", 3,
    "Primary 4", 4,
    "Primary 5", 5,
    "Primary 6", 6,
    "Primary 7", 7,
    "Primary 8", 8,
    "Secondary 1", 9,
    "Secondary 2", 10,
    "Secondary 3", 11,
    "Secondary 4", 12,
    "College", 13,
    "University", 13
  ) 
  cov_data = cov_data %>%
    mutate(
      completed_primary = (school == "Primary 8" | str_detect(school, "Secondary|College|University"))
    ) %>%
    left_join(school_year_df) %>%
    mutate(
      have_phone_lgl = case_when(
        have_phone == "Yes" ~ TRUE, 
        have_phone == "No" ~ FALSE, 
        TRUE ~ NA
      )
    )

  # ethnicity and religion
  cov_data = cov_data %>%
    mutate(
      ethnicity_luhya = ethnicity == "Luhya",
      religion_christianity = religion == 1
    )

    return(cov_data)
}


baseline_worm = baseline.data %>%
  clean_worm_covariates()
clean_takeup_variables = function(data) {
  data %>%
    mutate(
      female = gender == "female",
      have_phone_lgl = phone_owner
    )
}
analysis_data = analysis_data %>%
  clean_takeup_variables()

cluster_treat_df = read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))  %>%
  mutate(
      treat_dist = paste0(
      "treat: ", 
      assigned.treatment,
      ", dist: ", dist.pot.group
      ) %>% factor()
  ) %>%
  select(cluster.id, treat_dist, cluster.dist.to.pot = dist.to.own.pot) %>%
  unique()


pretreat_data = clean_pretreat_covariates(baseline.data, endline.data) %>%
  left_join(cluster_treat_df, by = "cluster.id") %>%
  filter(!is.na(treat_dist))

## Baseline Balance
baseline_worm_data = baseline_worm %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  filter(!is.na(treat_dist))

# monitoring checks
sens_imp_df = read_csv("data/raw-data/Sensitization Monitoring Form.csv")
sens_imp_hh_df = read_csv("data/raw-data/Sensitization Monitoring Form-household.csv")


unique_hh_message_df = sens_imp_hh_df %>%
  mutate(
    message_list = str_split(message, " ")
  ) %>%
  select(PARENT_KEY, message_list) %>%
  unnest(message_list) %>%
  group_by(PARENT_KEY) %>%
  unique() %>%
  summarise(message_list = list(as.numeric(message_list))) %>%
  mutate(
    knowledge_message = map_lgl(message_list, ~all(c(1, 2, 3, 4, 5) %in% .x)
     ),
    availability_message = map_lgl(message_list, ~all(c(6, 7) %in% .x))
    )



clean_sens_imp_df = sens_imp_df %>%
  filter(!is.na(enumerator)) %>%
  filter(!is.na(announcement)) %>%
  mutate(announce_church = str_detect(where, "1")) %>%
  left_join(unique_hh_message_df, by = c("KEY" = "PARENT_KEY")) %>%
  select(
    cluster.id = cluster_id,
    announcement,
    announce_church,
    knowledge_message,
    availability_message
  ) %>%
  inner_join(cluster_treat_df, by = "cluster.id") %>%
  filter(!is.na(treat_dist)) %>%
  left_join(
    cluster.strat.data %>%
      select(cluster.id, county)
  ) 


# PoT verification monitoring
clean_sens_summ_imp_df = clean_sens_imp_df %>%
  group_by(cluster.id, treat_dist, county) %>%
  summarise(
    n_announce = sum(announcement, na.rm = TRUE),
    n_announce_church = sum(announce_church, na.rm = TRUE),
    n_knowledge_message = sum(knowledge_message, na.rm = TRUE),
    n_avail_message = sum(availability_message, na.rm = TRUE),
    n_message = sum(!is.na(knowledge_message) | !is.na(availability_message)),
    n_total = n()
  ) %>%
  mutate(
    pct_announce = n_announce / n_total,
    pct_church = n_announce_church / n_total,
    pct_knowledge_message = n_knowledge_message / n_message,
    pct_avail_message = n_avail_message / n_message
  ) 
  
  


clean_implementation_vars = function(data)  {
  data %>%
    mutate(
      know_deworm = know_deworm == "yes",
      treat_begin = treat_begin == "knows",
      treat_end = treat_end == "knows",
      chv_visit =  chv_visit == "yes",
      flyer = flyer == "yes"
    )
}


worm_vars = c(
  "treated_lgl", 
  "know_how_stop_worms",
  "all_can_get_worms",
  "correct_when_treat",
  "fully_aware_externalities",
  "know_worms_infectious"
)


pretreat_vars = c(
  "completed_primary", 
  "floor_tile_cement",
  "ethnicity_luhya",
  "religion_christianity"
)

implementation_vars = c(
  "chv_visit",
  "flyer"
)


sens_vars = c(
  "pct_announce",
  "pct_church",
  "pct_knowledge_message",
  "pct_avail_message"
  )

# Distance/cluster size balance vars
takeup_vars = c(
  "n_per_cluster",
  "female",
  "have_phone_lgl",
  "age.census",
  "cluster.dist.to.pot"
)


# Indiv level balance variables
census_vars = c(
  "dewormed",
  "know_age" # just include this so fixest creates a list of fits...
)


praise_vars = c(
  "praise_immuniz",
  "praise_dewor"
)

stigma_vars = c(
  "stigma_immuniz",
  "stigma_dewor"
)


clean_census_data = census_data %>%
  filter(!is.na(assigned.treatment)) %>%
  right_join(
    analysis_data %>%
      select(cluster.id, dist.pot.group, cluster.dist.to.pot) %>%
      unique()
  ) %>%
  mutate(
    female = gender == 2,
    age = age.census,
    have_phone_lgl = have_phone == "Yes" 
    )

analysis_data = analysis_data %>%
    left_join(
        n_indiv_df,
        by = "cluster.id"
    ) %>%
    ungroup() %>%
    mutate(have_phone_lgl = have_phone == "Yes")



analysis_data = analysis_data %>%
  mutate(
      treat_dist = paste0(
      "treat: ", 
      assigned.treatment,
      ", dist: ", dist.pot.group
      ) %>% factor()
    )  


#### Functions for RF bootstrap ------------------------------------------------
# Estimate and Predict to generate ATEs
quick_pred = function(data) {
    fit = 
        feglm(
            dewormed ~ 0 + assigned_treatment +  standard_cluster.dist.to.pot   + i(assigned_treatment, standard_cluster.dist.to.pot, "control") | county,
            data = data, 
            family = binomial(link = "probit"),
            nthreads = 1
        )
    data$pred = predict(fit)
    signal_fit = 
        feglm(
            dewormed ~ 0 + signal +   standard_cluster.dist.to.pot + i(signal, standard_cluster.dist.to.pot, "no signal") | county,
            data = data, 
            family = binomial(link = "probit"),
            nthreads = 1
        )
    data$signal_pred = predict(signal_fit)
    data = data %>%
        select(
            assigned_dist_group,
            assigned_treatment,
            signal,
            standard_cluster.dist.to.pot,
            pred,
            signal_pred
        )
    return(data)
}

# For a given set of IDs, create bs data - n.b. this allows cluster to appear 
# multiple times
create_bs_data = function(split_data, ids) {
    bs_df = do.call(rbind, split_data[ids])
    return(bs_df)
}


# Create bootstrapped predictions
create_bs_preds = function(pred_df, ...) {
    preds = pred_df %>%
        group_by(
            assigned_treatment,
            assigned_dist_group,
            ...
        ) %>%
        summarise(
            mean_pred = mean(pred),
            .groups = "drop"
        ) %>%
        bind_rows(
            pred_df %>%
                group_by(
                    assigned_treatment,
                    ...
                ) %>%
                summarise(
                    mean_pred = mean(pred),
                    .groups = "drop"
                ) %>%
                mutate(
                    assigned_dist_group = "combined"
                )
        )
        signal_pred = pred_df %>%
            group_by(signal, assigned_dist_group, ...) %>%
            summarise(
                mean_pred = mean(signal_pred),
                .groups = "drop"
            ) %>%
            bind_rows(
            pred_df %>%
                group_by(
                    signal, 
                    ...
                ) %>%
                summarise(
                    mean_pred = mean(signal_pred)
                ) %>%
                mutate(
                    assigned_dist_group = "combined"
                )
            )
    preds = bind_rows(preds, signal_pred)
       return(preds)
}
bs_fit = function(seed, split_data = split_analysis_data, ...) {
    set.seed(seed)
    ids = names(split_data)
    sampled_ids = sample(ids, length(ids), replace = TRUE)
    bs_df = create_bs_data(split_data, sampled_ids)
    bs_fit = quick_pred(bs_df) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
}


# Function to generate samples from a Dirichlet distribution
generate_dirichlet <- function(alpha, n) {
  # Generate matrix of gamma samples
  gamma_samples <- matrix(rgamma(n * length(alpha), shape = alpha, scale = 1), nrow = n)
  
  # Normalize rows to sum to 1
  dirichlet_samples <- gamma_samples / rowSums(gamma_samples)
  
  return(dirichlet_samples)
}

pred_bs_f = function(f, f_signal, data, weights, realised_fit = FALSE) {
    if (realised_fit == TRUE) {
        data$wt = 1
    } else {
        data$wt = weights[data$cluster_id]
    }
    fit = f(data, weights = ~wt)
    data$pred = predict(fit)
    signal_fit = f_signal(data, weights = ~wt)
    data$signal_pred = predict(signal_fit)
    data = data %>%
        select(
            assigned_dist_group,
            assigned_treatment,
            signal,
            standard_cluster.dist.to.pot,
            pred,
            signal_pred,
            any_of("sms_treatment")
        )
    return(data)
}





bayes_bs_f = function(seed, f, f_signal, data, ...) {
    set.seed(seed)
    n_clusters = length(unique(data$cluster.id))
    alpha = rep(1, n_clusters)
    weights = generate_dirichlet(alpha, 1)
    bs_fit = pred_bs_f(f, f_signal, data, weights = weights) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 

actual_bayesian_bs_fit = function(seed, f, f_signal, data, ...) {
    bs_fit = pred_bs_f(f, f_signal, data, 1, realised_fit = TRUE) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 


add_predictions = function(draws, ...) {
    bra_minus_cal = draws %>%
            filter(assigned_treatment == "bracelet" | assigned_treatment == "calendar") %>%
            group_by(seed, assigned_dist_group, ...) %>%
            summarise(
                mean_pred = mean_pred[assigned_treatment == "bracelet"] - mean_pred[assigned_treatment == "calendar"],
            ) %>%
            mutate(
                assigned_treatment = "bracelet - calendar"
            )
    
    bra_minus_control = draws %>%
            filter(
              assigned_treatment == "bracelet" 
              ) %>%
            group_by(seed, assigned_dist_group, ...) %>%
            summarise(
                mean_pred = mean_pred[assigned_treatment == "bracelet"],
            ) %>%
            mutate(
                assigned_treatment = "bracelet - control"
            )

    bra_minus_no_signal = bind_rows(
      bra_minus_cal,
      bra_minus_control
    ) %>%
    group_by(seed, assigned_dist_group, ...) %>%
    summarise(
      mean_pred = mean(mean_pred)
    ) %>%
    mutate(assigned_treatment = "bracelet - no signal") 

    draws = bind_rows(draws, bra_minus_cal, bra_minus_no_signal)

    fc_draws = draws %>%
    group_by(
        seed,
        assigned_treatment, ...
    ) %>%
    summarise(
        mean_pred = mean_pred[assigned_dist_group == "far"] - mean_pred[assigned_dist_group == "close"],
        assigned_dist_group = "far - close",
        .groups = "drop"
    )

    draws = bind_rows(draws, fc_draws) %>%
    ungroup() %>%
    mutate(
        assigned_treatment = factor(
            assigned_treatment,
            levels = c("control", "bracelet - calendar", "bracelet - no signal", "ink", "calendar", "bracelet")
        ),
        assigned_dist_group = factor(
            assigned_dist_group,
            levels = c("combined", "close", "far", "far - close")
        )
    ) 

    return(draws)
}

add_signal_predictions = function(draws, ...) {
    fc_draws = draws %>%
    group_by(
        seed,
        signal,
        ...
    ) %>%
    summarise(
        mean_pred = mean_pred[assigned_dist_group == "far"] - mean_pred[assigned_dist_group == "close"],
        assigned_dist_group = "far - close"
    )
    draws = bind_rows(draws, fc_draws) %>%
        ungroup() %>%
        mutate(
            signal = factor(
                signal,
                levels = c("no signal", "signal")
            ),
            assigned_dist_group = factor(
                assigned_dist_group,
                levels = c("combined", "close", "far", "far - close")
            )
        ) 


    return(draws)
}

create_tes = function(draws, ...) {
    draws %>%
        group_by(seed, assigned_dist_group, ...) %>%
        mutate(
            mean_pred = if_else(assigned_treatment == "control", mean_pred, mean_pred - mean_pred[assigned_treatment == "control"])
        )
}

create_signal_tes = function(draws, ...) {
    draws %>%
        group_by(seed, assigned_dist_group, ...) %>%
        mutate(
            mean_pred = if_else(signal == "no signal", mean_pred, mean_pred - mean_pred[signal == "no signal"])
        )
}


clean_signal_draws = function(draws, ...) {
    draws %>%
        filter(!is.na(signal)) %>%
        select(-assigned_treatment) %>%
        create_signal_tes(., ...) %>%
        add_signal_predictions(., ...) %>%
        rename(
            assigned_treatment = signal,
            estimate = mean_pred
        )
}

clean_te_draws = function(draws, ...) {
    draws %>%
        filter(!is.na(assigned_treatment)) %>%
        select(-signal) %>%
        create_tes(., ...) %>%
        add_predictions(., ...)  %>%
        rename(estimate = mean_pred)
}

estimate_actual_fit = function(split_data = split_analysis_data) {
    ids = names(split_data)
    sampled_ids = ids
    bs_df = create_bs_data(split_data, sampled_ids)
    bs_fit = quick_pred(bs_df) %>%
        create_bs_preds()
    bs_fit$seed = "actual fit"
    return(bs_fit)
}


round_pval = function(pvals, digits = 3) {
    pvals = round(pvals, digits)
    pvals = if_else(pvals == 0, "<0.001", as.character(pvals))
    return(pvals)
}

add_summ_stats = function(bs_draws, actual_fit, ci_width = 0.95) {
    clean_tes = bs_draws %>%
      group_by(
          assigned_treatment,
          assigned_dist_group,
      ) %>%
      summarise(
          std_error = sd(estimate),
          conf.low = quantile(estimate, (1 - ci_width)/2),
          conf.high = quantile(estimate, 1 - (1 - ci_width)/2)
      ) %>%
      left_join(
          actual_fit,
          by = c("assigned_dist_group", "assigned_treatment")
      ) %>%
      mutate(
          pval = 2*pnorm(-abs(realised_pred)/std_error),
          oneside_pval = pnorm(-realised_pred/std_error)
      ) %>%
      mutate(
          pval = round_pval(pval, 3),
          oneside_pval = round_pval(oneside_pval, 3)
      ) %>%
      select(
          assigned_treatment, 
          assigned_dist_group, 
          realised_pred, 
          std_error, 
          conf.low,
          conf.high,
          pval, 
          oneside_pval) %>%
      rename(estimate = realised_pred) 
      return(clean_tes)
}

# wrapper function for all of the above
create_regression_output = function(data, f, f_signal, B_draws = 500, 
                                    stat = params$stat,
                                    caption = "Average Treatment Effects: Reduced Form",
                                    dependent_var = "Dependent variable: Take-up") {
  bs_draws = map_dfr(
    1:B_draws,
    ~bayes_bs_f(
      seed = .x,
      f = f,
      f_signal = f_signal,
      data = data
    ),
    .progress = TRUE
    )
  
  clean_te_draws = bs_draws %>%
    clean_te_draws()
  clean_signal_draws = bs_draws %>%
    clean_signal_draws()

  realised_fit = actual_bayesian_bs_fit(
    seed = "realised fit",
    f = f,
    f_signal = f_signal,
    data = data
  )
  signal_fit = realised_fit %>%
    clean_signal_draws() %>%
    rename(realised_pred = estimate) %>%
    select(realised_pred, assigned_dist_group, assigned_treatment)
  te_fit = realised_fit %>%
    clean_te_draws() %>%
    rename(realised_pred = estimate) %>%
    select(realised_pred, assigned_dist_group, assigned_treatment)

  signal_summ = add_summ_stats(clean_signal_draws, signal_fit)
  te_summ = add_summ_stats(clean_te_draws, te_fit)


  pval_only_terms = c("bracelet - calendar", "signal")

  overall_summ = bind_rows(
    signal_summ,
    te_summ
  ) %>%
    mutate(
      show_pval_only = assigned_treatment %in% pval_only_terms
    ) %>%
    filter(assigned_treatment != "no signal") 


  default_tbl = overall_summ %>%
    prep_tbl(stat = params$stat) %>%
    nice_kbl_table(
      cap = caption,
      outcome_var = dependent_var
      )


  different_order_tbl = overall_summ %>%
    prep_tbl(stat = stat) %>%
    mutate(
      assigned_treatment = fct_relevel(
        assigned_treatment, 
        c(
          "Control", 
          "Bracelet - No Signal", 
          "$H0$: Any Signal > No Signal, $p$-value",
          "$H0$: Any Signal $\\neq$ No Signal, $p$-value",
          "$H0$: Bracelet > Calendar, $p$-value",
          "$H0$: Bracelet $\\neq$ Calendar, $p$-value",
          "Bracelet", "Calendar", "Ink"
        ))) %>% 
      arrange(assigned_treatment) %>%
    nice_kbl_table(
      cap = caption,
      outcome_var = dependent_var
    )

  return(list(
    tidy_summary = overall_summ,
    default_tbl = default_tbl,
    different_order_tbl = different_order_tbl
  ))

}



### Table/Kable functions ------------------------------------------------------
prep_tbl = function(tes, stat = "ci") {

    tbl_dist_levels = c(
        "combined",
        "close",
        "far",
        "far - close"
    )

    tbl_contrast_levels = c(
        "bracelet",
        "calendar",
        "ink",
        "control",
        "bracelet - no signal",
        "signal",
        "bracelet - calendar",
        "signal two-side pval",
        "bracelet - calendar two-side pval"
    )



    if (stat == "ci") {
        tes = tes %>%
            mutate(
                val = paste0(
                    "(",
                    round(conf.low, 3),
                    ", ",
                    round(conf.high, 3),
                    ")"
                )
            )
    } else if (stat == "std.error"){
        tes = tes %>%
            mutate(val = paste0("{[", round(std_error, 3), "]}"))
    } else {
        tes = tes %>%
            mutate(
                val = paste0("{[", round_pval(pval, 3), "]}")
            )
    }





    tbl =  tes %>%
        select(assigned_treatment, assigned_dist_group, estimate, conf.low, conf.high, val, pval, oneside_pval, show_pval_only)  %>%
        mutate(across(where(is.numeric), ~round(.x, 3))) %>%
        mutate(
            estim_std = if_else(
                show_pval_only,
                linebreak(paste0(oneside_pval), align = "c"),
                linebreak(paste0(estimate,"\n", str_glue("{val}")), align = "c") 
            )
        ) %>%
        bind_rows(
            filter(., show_pval_only == TRUE) %>%
              mutate(
                estim_std = linebreak(paste0(pval), align = "c")
              ) %>%
              mutate(
                assigned_treatment = paste0(assigned_treatment, " two-side pval")
              )
        ) %>%
        select(assigned_treatment, assigned_dist_group, estim_std) %>%
        mutate(
            assigned_dist_group = factor(assigned_dist_group, tbl_dist_levels),
            assigned_dist_group = fct_relabel(assigned_dist_group, str_to_title),
            assigned_treatment = factor(assigned_treatment, tbl_contrast_levels)
        ) %>%
        arrange(assigned_dist_group, assigned_treatment) %>%
        pivot_wider(
            names_from = assigned_dist_group,
            values_from = estim_std
        ) %>%
        mutate(
            assigned_treatment = fct_relabel(assigned_treatment, str_to_title),
            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Any Signal > No Signal, $p$-value"  = "Signal"), 
            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Bracelet > Calendar, $p$-value" = "Bracelet - Calendar"),


            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Any Signal $\\neq$ No Signal, $p$-value"  = "Signal Two-Side Pval"), 
            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Bracelet $\\neq$ Calendar, $p$-value" = "Bracelet - Calendar Two-Side Pval")
        ) 
    
    return(tbl)
}

nice_kbl_table = function(tbl, cap, outcome_var = "Dependent variable: Take-up", stat = params$stat) {
  linesep_str = if_else(stat == "ci", "\\addlinespace", "")

  nice_kbl = tbl %>%
  kbl(
    col.names = c(
      # "Estimand", 
      # "Treatment", 
      outcome_var,
      paste0("(", 1:4, ")")
    ), 
    format = "latex", 
    linesep = linesep_str, 
    booktabs = TRUE, 
    escape = FALSE, 
    align = "lcccc", 
    caption = cap
  )  %>%
  kable_styling(
    latex_options = c("scale_down")
  ) %>%
  add_header_above(
    c(" ", 
      "Combined", 
      "Close", 
      "Far", 
      "Far - Close"
      ), 
    line = FALSE
  ) %>%
  add_header_above(
    c(
      " " = 1,
      "Reduced Form" = 4
      )
  ) %>%
  row_spec(c(6), hline_after = TRUE) 
}

custom_save_latex_table = function(table, table_name, table_output_path = params$table_output_path){
  table_conn = file(
    file.path(
      table_output_path, paste0(table_name, ".tex")
    )
  )
  # Create space in latex due to rmd bug with \\ immediately followed by [
  # table = table  %>%
  #   str_replace_all(., "removeme12345", "  ")

  attr(table, "kable_meta")$contents = str_replace_all(attr(table, "kable_meta")$contents, "removeme12345", " ")
  table[1] = str_replace_all(table[1], "removeme12345", " ")


  clean_table = table %>%
    str_remove(
      ., 
      fixed("\\begin{table}")
    ) %>%
    str_remove(
      .,
      "\\\\caption\\{.*\\}"
    ) %>%
    str_remove(
      ., 
      "\\\\end\\{table\\}"
    ) 
    
    
    clean_table %>%
      writeLines(
        table_conn
      )
    close(table_conn)

    return(table)
}

#### Frequentist Estimates ####
split_analysis_data = split(analysis_data, analysis_data$cluster.id)
#| freq-estimates
analysis_data = analysis_data %>%
    mutate(
        county = factor(county),
        cluster.id = factor(cluster.id),
        assigned_treatment = assigned.treatment,
        assigned_dist_group = dist.pot.group,
        signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
        signal = factor(signal, levels = c("no signal", "signal"))
    )


cluster_dispersion_df = analysis_data %>%
  group_by(
    assigned_treatment,
    cluster_id
  ) %>%
  summarise(
    mse_dist_to_cluster = mean(((dist.to.pot - cluster.dist.to.pot)/1000)^2)
  ) %>%
  mutate(
    dispersed_community = mse_dist_to_cluster > 0.5
  ) %>%
  ungroup()

outlier_analysis_data = analysis_data %>%
    left_join(cluster_dispersion_df %>% select(-assigned_treatment), by = "cluster_id")  

no_outlier_analysis_data = outlier_analysis_data %>%
  filter(!dispersed_community) %>%
  group_by(cluster.id) %>%
  mutate(cluster_id = cur_group_id()) %>%
  ungroup()

# not super kosher but additional robustness check we can perform if refs ask
# splitting some dispersed clusters into a close/far cluster
# split_outlier_data = outlier_analysis_data %>%
#   group_by(
#     cluster_id
#   ) %>%
#   mutate(
#     mean_dist.to.pot = mean(dist.to.pot),
#     cluster_split = case_when(
#       dispersed_community == FALSE ~ as.character(cluster_id), 
#       dispersed_community == TRUE & dist.to.pot < mean_dist.to.pot ~ paste0(cluster_id, " - close"),
#       dispersed_community == TRUE & dist.to.pot >= mean_dist.to.pot ~ paste0(cluster_id, " - far")
#     )
#   )  %>%
#   ungroup() %>%
#   group_by(cluster_split) %>%
#   mutate(mean_dist.to.pot = mean(dist.to.pot)) %>%
#   mutate(
#     cluster_id = cur_group_id(),
#     assigned_dist_group = if_else(str_detect(cluster_split, "close") & mean_dist.to.pot < 1250, "close", assigned_dist_group),
#     assigned_dist_group = if_else(str_detect(cluster_split, "far") & mean_dist.to.pot > 1250, "far", assigned_dist_group)
#   )  %>%
#   ungroup()

#### Distance Checks -----------------------------------------------------------
outlier_analysis_data %>%
  filter(dispersed_community == TRUE) %>%
  select(
    cluster_id,
    assigned_treatment,
    cluster.dist.to.pot,
    dist.to.pot
  ) %>%
  group_by(cluster_id) %>%
  mutate(
    mean_dist.to.pot = mean(dist.to.pot),
    median_dist.to.pot = median(dist.to.pot)
  ) %>%
  ggplot() +
  geom_histogram(
    aes(x = dist.to.pot, fill = factor(cluster_id))
  ) +
  geom_vline(
    aes(
      xintercept = mean_dist.to.pot
    ),
    linetype = "dotted"
  ) +
  geom_vline(
    aes(
      xintercept = cluster.dist.to.pot
    ),
    linetype = "longdash"
  ) +
  geom_vline(
    aes(
      xintercept = median_dist.to.pot
    ),
    linetype = "dotdash"
  ) +
  facet_wrap(~assigned_treatment + cluster_id) +
  labs(
    x = "Distance to PoT (m)",
    y = "Count",
    fill = "Cluster",
    title = "Distribution of Distances in Outlier Clusters",
    caption = "Dashed line represents the cluster's centroid distance to the PoT. Dotted line the mean distance to the PoT."
  )
# ggsave(
#   "temp-data/dist-distance-to-pot-outliers.png",
#   width = 10,
#   height = 10
# )
#### Distance Checks End -------------------------------------------------------
main_spec_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control")  | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}

main_spec_signal_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + signal + standard_cluster.dist.to.pot + i(signal, standard_cluster.dist.to.pot, "no signal") | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}



main_spec_output = create_regression_output(
  data = analysis_data,
  f = main_spec_regression,
  f_signal = main_spec_signal_regression
)

main_spec_output$tidy_summary %>%
  write_csv("temp-data/reducedform-tidy-tes.csv")  


main_spec_output$default_tbl %>%
  custom_save_latex_table(
    table_name = "rf_main_spec_tbl"
  )

main_spec_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_main_spec_tbl_weird_order"
  )


# coefs on main specification interaction terms
main_fit = feglm(
    dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control")  | county, 
    data = analysis_data,
    family = binomial(link = "probit"),
    nthreads = 1,
    cluster = ~cluster.id
)

etable(
  main_fit,
  dict = c(
    "assigned_treatment = ink" = "Ink",
    "assigned_treatment = calendar" = "Calendar",
    "assigned_treatment = bracelet" = "Bracelet",
    assigned_treatmentink = "Ink",
    assigned_treatmentcalendar = "Calendar",
    assigned_treatmentbracelet = "Bracelet",
    "standard_cluster.dist.to.pot" = "Distance to PoT",
    "county" = "County",
    "cluster.id" = "Cluster",
    dewormed = "Dewormed"
    ),
    headers = "",
    replace = TRUE
  )


etable(
  main_fit,
  dict = c(
    "assigned\\_treatment $=$ ink" = "Ink",
    "assigned\\_treatment = calendar" = "Calendar",
    "assigned\\_treatment = bracelet" = "Bracelet",
    assigned_treatmentink = "Ink",
    assigned_treatmentcalendar = "Calendar",
    assigned_treatmentbracelet = "Bracelet",
    "standard_cluster.dist.to.pot" = "Distance to PoT",
    "county" = "County",
    "cluster.id" = "Cluster",
    dewormed = "Dewormed",
    "assigned_treatment" = "Treatment",
    "ink" = "Ink",
    "calendar" = "Calendar",
    "bracelet" = "Bracelet"
    ),
    headers = "",
    tex = TRUE,
    style.tex = style.tex("aer"),
    replace = TRUE,
    file = 
    file.path(
      params$table_output_path, paste0("main_spec_regression_coefs", ".tex")
    )
  )



# main specification levels
main_spec_bs_draws = map_dfr(
  1:500,
  ~bayes_bs_f(
    seed = .x,
    f = main_spec_regression,
    f_signal = main_spec_signal_regression,
    data = analysis_data
  ),
  .progress = TRUE
  )

main_spec_levels = actual_bayesian_bs_fit(
  seed = "realised fit",
  f = main_spec_regression,
  f_signal = main_spec_signal_regression,
  data = analysis_data
) %>%
  filter(!is.na(assigned_treatment)) 

main_spec_levels_ci = main_spec_bs_draws %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  summarise(
    conf.low = quantile(mean_pred, 0.025),
    conf.high = quantile(mean_pred, 0.975)
  ) %>%
  filter(!is.na(assigned_treatment))

tidy_main_spec_levels = left_join(
  main_spec_levels,
  main_spec_levels_ci,
  by = c("assigned_treatment", "assigned_dist_group")
) %>%
  select(-signal, -seed) %>%
  rename(estimate = mean_pred)
tidy_main_spec_levels %>%
  write_csv("temp-data/reducedform-tidy-levels.csv")  

# Distance entering with its square
nonlinear_distance_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + standard_cluster.dist.to.pot^2 + i(assigned_treatment, standard_cluster.dist.to.pot, "control") + i(assigned_treatment, standard_cluster.dist.to.pot^2, "control") | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}

nonlinear_distance_signal_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + signal + standard_cluster.dist.to.pot + standard_cluster.dist.to.pot^2 + i(signal, standard_cluster.dist.to.pot, "no signal") + i(signal, standard_cluster.dist.to.pot^2, "no signal") | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}
nonlinear_distance_output = create_regression_output(
  data = analysis_data,
  f = nonlinear_distance_regression,
  f_signal = nonlinear_distance_signal_regression
)

nonlinear_distance_output$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-nonlinear-dist-tidy-tes.csv")  


nonlinear_distance_output$default_tbl %>%
  custom_save_latex_table(
    table_name = "rf_nonlinear_dist_tbl"
  )

nonlinear_distance_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_nonlinear_dist_tbl_weird_order"
  )

# discrete distance
discrete_distance_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + assigned_treatment + assigned_dist_group  + i(assigned_treatment, assigned_dist_group, "control")  | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}
discrete_distance_signal_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + signal + assigned_dist_group  + i(signal, assigned_dist_group, "no signal")  | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}

discrete_distance_output = create_regression_output(
  data = analysis_data,
  f = discrete_distance_regression,
  f_signal = discrete_distance_signal_regression
)

discrete_distance_output$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-discrete-dist-tidy-tes.csv")  


discrete_distance_output$default_tbl %>%
  custom_save_latex_table(
    table_name = "rf_discrete_dist_tbl"
  )

discrete_distance_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_discrete_dist_tbl_weird_order"
  )

## HH Dist regression
hh_spec_regression = function(data, weights) {
  feglm(
    dewormed ~  0  + assigned_treatment + dist.to.pot + i(assigned_treatment, dist.to.pot, "control") | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}
hh_spec_signal_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + signal + dist.to.pot + i(signal, dist.to.pot, "no signal") | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}

hh_spec_output = create_regression_output(
  data = analysis_data,
  f = hh_spec_regression,
  f_signal = hh_spec_signal_regression
)


hh_spec_output$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-hhdist-tidy-tes.csv")  


hh_spec_output$default_tbl %>%
  custom_save_latex_table(
    table_name = "rf_hh_spec_tbl"
  )

hh_spec_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_hh_spec_tbl_weird_order"
  )


# community dist + controlling for HH dist
community_control_spec_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") + dist.to.pot | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}
community_control_spec_signal_regression = function(data, weights) {
  feglm(
    dewormed ~ 0 + signal + standard_cluster.dist.to.pot  + i(signal, standard_cluster.dist.to.pot, "no signal") + dist.to.pot | county, 
    data = data,
    family = binomial(link = "probit"),
    nthreads = 1,
    weights = ~wt
  )
}

community_control_spec_output = create_regression_output(
  data = analysis_data,
  f = community_control_spec_regression,
  f_signal = community_control_spec_signal_regression
)


community_control_spec_output$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-communitycontrol-tidy-tes.csv")  



community_control_spec_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_communitycontrol_spec_tbl_weird_order"
  )

#### Outliers-------------------------------------------------------------------

no_outlier_spec = create_regression_output(
  data = no_outlier_analysis_data,
  f = main_spec_regression,
  f_signal = main_spec_signal_regression
)

no_outlier_spec$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-nooutlier-tidy-tes.csv")  


no_outlier_spec$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_nooutlier_spec_tbl_weird_order"
  )

## Community dist + HH dist control with no outliers
no_outlier_community_control_spec_output = create_regression_output(
  data = no_outlier_analysis_data,
  f = community_control_spec_regression,
  f_signal = community_control_spec_signal_regression
)
no_outlier_community_control_spec_output$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-nooutliercommunitycontrol-tidy-tes.csv")  
no_outlier_community_control_spec_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_nooutliercommunitycontrol_spec_tbl_weird_order"
  )
#### Beliefs -------------------------------------------------------------------

disagg_base_belief_data = analysis_data %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  nest_join(
    endline.know.table.data %>% 
      filter(fct_match(know.table.type, "table.A")),
    by = "KEY.individ", 
    name = "knowledge_data"
  ) %>% 
  mutate(
    map_dfr(knowledge_data, ~ {
      tibble(
        obs_know_person = sum(.x$num.recognized),
        obs_know_person_prop = mean(.x$num.recognized),
        knows_other_dewormed = sum(fct_match(.x$dewormed, c("yes", "no")), na.rm = TRUE),
        knows_other_dewormed_yes = sum(fct_match(.x$dewormed, "yes"), na.rm = TRUE),
        knows_other_dewormed_no = sum(fct_match(.x$dewormed, "no"), na.rm = TRUE),
        thinks_other_knows = sum(fct_match(.x$second.order, c("yes", "no")), na.rm = TRUE),
        thinks_other_knows_yes = sum(fct_match(.x$second.order, "yes"), na.rm = TRUE),
        thinks_other_knows_no = sum(fct_match(.x$second.order, "no"), na.rm = TRUE),
      )
    }
  )) %>%
    filter(obs_know_person > 0)  %>%
    select(
      KEY.individ, 
      contains("know"), 
      assigned.treatment, 
      dist.pot.group, 
      assigned_dist_group,
      cluster.id,
      cluster.dist.to.pot,
      standard_cluster.dist.to.pot,
      dist.to.pot,
      county
      ) %>%
    mutate(
        doesnt_know_other_dewormed = obs_know_person - knows_other_dewormed, 
        doesnt_think_other_knows = obs_know_person - thinks_other_knows
    ) %>% 
    select(KEY.individ, 
           assigned.treatment,
           assigned_dist_group,
           obs_know_person,
           knows_other_dewormed_yes,
           knows_other_dewormed_no,
           doesnt_know_other_dewormed, 
           thinks_other_knows_yes, 
           thinks_other_knows_no, 
           doesnt_think_other_knows,
           cluster.id,
           cluster.dist.to.pot,
           standard_cluster.dist.to.pot,
           dist.to.pot,
           county
           ) %>%
    gather(variable, value, 
        knows_other_dewormed_yes:doesnt_think_other_knows)   %>%
    mutate(knowledge_type = case_when(
        str_detect(variable, "_yes") ~ "yes",
        str_detect(variable, "_no") ~ "no",
        str_detect(variable, "doesnt") ~ "doesn't know"
    )) %>%
    mutate(belief_type = if_else(str_detect(variable, "think"), "2ord", "1ord")) %>%
    mutate(prop = value/obs_know_person) 

know_df = disagg_base_belief_data %>%
  filter(knowledge_type == "doesn't know") %>%
  mutate(
    prop_knows = 1 - prop
  ) %>%
  group_by(cluster.id) %>%
  mutate(cluster_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(
      county = factor(county),
      cluster.id = factor(cluster.id),
      assigned_treatment = assigned.treatment,
      signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
      signal = factor(signal, levels = c("no signal", "signal"))
  ) 
  




f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") | county,
    data = data,
    weights = weights
  )
}

f_know_signal = function(data, weights) {
  feols(
    prop_knows ~ signal + standard_cluster.dist.to.pot + i(signal, standard_cluster.dist.to.pot, "no signal") | county,
    data = data,
    weights = weights
  )
}

hh_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + dist.to.pot + i(assigned_treatment, dist.to.pot, "control") | county,
    data = data,
    weights = weights
  )
}

hh_f_know_signal = function(data, weights) {
  feols(
    prop_knows ~ signal + dist.to.pot + i(signal, dist.to.pot, "no signal") | county,
    data = data,
    weights = weights
  )
}


discrete_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") | county,
    data = data,
    weights = weights
  )
}
discrete_f_know_signal = function(data, weights) {
  feols(
    prop_knows ~ signal + assigned_dist_group + i(signal, assigned_dist_group, "no signal") | county,
    data = data,
    weights = weights
  )
}

discrete_fob_output = create_regression_output(
  data = know_df %>%
    filter(belief_type == "1ord"),
  f = discrete_f_know,
  f_signal = discrete_f_know_signal,
  dependent_var = "Dependent variable: First-order beliefs"
)
discrete_fob_output$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-discrete-fob-tidy-tes.csv")  
discrete_fob_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_discrete_fob_spec_tbl_weird_order"
  )

## SOB Main Spec ---------------------------------------------------------------
sob_fit = create_regression_output(
  data = know_df %>%
    filter(belief_type == "2ord"),
  f = f_know,
  f_signal = f_know_signal,
  dependent_var = "Dependent variable: Second-order beliefs"
)

sob_fit$tidy_summary %>%
  write_csv("temp-data/reducedform-sob-tidy-tes.csv")  
sob_fit$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_sob_spec_tbl_weird_order"
  )


## robustness HH dist
robust_hh_fob_output = create_regression_output(
  data = know_df %>%
    filter(belief_type == "1ord"),
  f = hh_f_know,
  f_signal = hh_f_know_signal,
  dependent_var = "Dependent variable: First-order beliefs"
)
robust_hh_fob_output$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-hhdist-fob-tidy-tes.csv")  
robust_hh_fob_output$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_hhdist_fob_spec_tbl_weird_order"
  )
 


fob_know_bs_draws = map_dfr(
    1:500,
    ~bayes_bs_f(
        seed = .x, 
        f = f_know, 
        f_signal = f_know_signal, 
        data = know_df %>% 
          filter(belief_type == "1ord")
    ),
    .progress = TRUE
)


know_bs_te_draws = fob_know_bs_draws %>%
  clean_te_draws()

know_bs_signal_draws = fob_know_bs_draws %>%
  clean_signal_draws()



realised_know_fit = actual_bayesian_bs_fit(
  seed = "realised fit", 
  f = f_know, 
  f_signal = f_know_signal,
  data = know_df %>% 
    filter(belief_type == "1ord"))

realised_know_signal_fit = realised_know_fit %>%
    clean_signal_draws() %>%
    rename(realised_pred = estimate) %>%
    select(assigned_dist_group, assigned_treatment, realised_pred)


realised_know_te_fit = realised_know_fit %>%
    clean_te_draws() %>%
    rename(realised_pred = estimate) %>%
    select(assigned_dist_group, assigned_treatment, realised_pred)


clean_know_signal_tes = add_summ_stats(know_bs_signal_draws, realised_know_signal_fit)
clean_know_tes = add_summ_stats(know_bs_te_draws, realised_know_te_fit)

pval_only_terms = c("bracelet - calendar", "signal")
clean_know_df = bind_rows(
  clean_know_tes,
  clean_know_signal_tes
) %>%
  mutate(
    show_pval_only = assigned_treatment %in% pval_only_terms
  ) %>%
  filter(assigned_treatment != "no signal")


clean_know_df %>%
  prep_tbl(stat = params$stat) %>%
  nice_kbl_table(
    cap = "Average Treatment Effects: Knowledge",
    outcome_var = "Dependent variable: First-order beliefs"
  ) %>%
  custom_save_latex_table(
    table_name = "rf_know_spec_tbl"
  )

clean_know_df %>%
  prep_tbl(stat = params$stat) %>%
  mutate(
    assigned_treatment = fct_relevel(
      assigned_treatment, 
      c(

        "Control", 
        "Bracelet - No Signal",
        "$H0$: Any Signal > No Signal, $p$-value",
        "$H0$: Any Signal $\\neq$ No Signal, $p$-value",
        "$H0$: Bracelet > Calendar, $p$-value",
        "$H0$: Bracelet $\\neq$ Calendar, $p$-value",
        "Bracelet", "Calendar", "Ink"
      ))) %>% 
    arrange(assigned_treatment) %>%
  nice_kbl_table(
    cap = "Average Treatment Effects: Knowledge",
    outcome_var = "Dependent variable: First-order beliefs"
  ) %>%
  custom_save_latex_table(
    table_name = "rf_know_spec_tbl_weird_order"
  )


    

clean_know_df %>%
  write_csv("temp-data/knowledge-tidy-tes.csv")  

#### FOB Main Spec + Controlling for HH Dist
f_know_control_hh = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + standard_cluster.dist.to.pot + dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") | county,
    data = data,
    weights = weights
  )
}

f_know_control_hh_signal = function(data, weights) {
  feols(
    prop_knows ~ signal + standard_cluster.dist.to.pot + dist.to.pot + i(signal, standard_cluster.dist.to.pot, "no signal") | county,
    data = data,
    weights = weights
  )
}


fob_control_hh = create_regression_output(
  data = know_df %>%
    filter(belief_type == "1ord"),
  f = f_know_control_hh,
  f_signal = f_know_control_hh_signal
)


fob_control_hh$tidy_summary %>%
  write_csv("temp-data/reducedform-robustness-fob-controlhh-tidy-tes.csv")  



fob_control_hh$different_order_tbl %>%
  custom_save_latex_table(
    table_name = "rf_fob_controlhh_spec_tbl_weird_order"
  )

## FOB Levels
# main specification levels
fob_bs_draws = map_dfr(
  1:500,
  ~bayes_bs_f(
    seed = .x,
    f = f_know,
    f_signal = f_know_signal,
    data = know_df %>%
      filter(belief_type == "1ord")
  ),
  .progress = TRUE
  )

fob_levels_point = actual_bayesian_bs_fit(
  seed = "realised fit",
  f = f_know,
  f_signal = f_know_signal,
  data = know_df %>%
    filter(belief_type == "1ord")
) %>%
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


#### SMS -----------------------------------------------------------------------


monitored_sms_data <- analysis.data %>% 
  filter(mon_status == "monitored") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()





sms_analysis_data <- monitored_sms_data %>% 
    mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group, 
    sms_treatment = sms.treatment.2, 
    phone_owner = if_else(phone_owner == TRUE, "phone", "nophone"), 
    sms_treatment = str_replace_all(sms_treatment, "\\.", "")) %>%
    # reminder.only only present in control condition
    filter(phone_owner == "phone") %>%
    mutate(sms_treatment = factor(sms_treatment)) %>%
    mutate(
        county = factor(county),
        cluster.id = factor(cluster.id),
        assigned_treatment = assigned.treatment,
        assigned_dist_group = dist.pot.group,
        signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
        signal = factor(signal, levels = c("no signal", "signal"))
    )





f_sms = function(data, weights) {
  feglm(
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
    weights = weights,
    family = binomial(link = "probit")
  )
}

f_sms_signal = function(data, weights) {
  feglm(
    dewormed ~ 0 +
      signal + 
      standard_cluster.dist.to.pot + 
      sms_treatment + 
      i(signal, standard_cluster.dist.to.pot, "no signal") +
      i(signal, sms_treatment, "no signal") +
      i(sms_treatment, standard_cluster.dist.to.pot) +
      sms_treatment:signal:standard_cluster.dist.to.pot 
      | county,
    data = data,
    weights = weights,
    family = binomial(link = "probit")
  )
}


sms_bs_draws = map_dfr(
    1:500,
    ~bayes_bs_f(
        seed = .x, 
        f = f_sms, 
        f_signal = f_sms_signal, 
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
  f_signal = f_sms_signal,
  data = sms_analysis_data,
  sms_treatment
)


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


clean_sms_tes %>%
  filter(sms_treatment != "smscontrol")  %>%
  mutate(show_pval_only = FALSE) %>%
  filter(sms_treatment != "reminderonly") %>%
  mutate(
    show_pval_only = assigned_treatment %in% pval_only_terms
  ) %>%
  prep_tbl(stat = params$stat) %>%
  nice_kbl_table(
    cap = "Heterogeneous SMS Average Treatment Effects",
    outcome_var = "Dependent variable: Take-up"
  ) %>%
  custom_save_latex_table(
    table_name = "sms_diff_tes_tbl"
  )




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

ggsave("temp-data/p-sms-tes.pdf", width = 8, height = 6)
