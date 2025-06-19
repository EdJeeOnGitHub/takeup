library(tidyverse)
library(broom)
library(data.table)
library(kableExtra)
library(knitr)
library(fixest)

if (interactive()) {
    params = lst(
        table_output_path = "presentations/rf-tables/main-specs",
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


dir.create("presentations/rf-tables/main-specs", showWarnings = FALSE)


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


## Loading Scripts
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"
## Loading Data
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

analysis_data %>%
  write_csv(
    "temp-data/analysis-data.csv"
  )

sd_of_dist = sd(analysis_data$cluster.dist.to.pot)
dist_sd = sd(analysis_data$cluster.dist.to.pot)




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
## Cleaning up analysis data
analysis_data = analysis_data %>%
  clean_takeup_variables()
## Getting cluster treatment assignment
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


colnames(endline.data)


endline.data %>%
  group_by(dist.pot.group) %>%
  summarise(
    mean_time_travel = mean(time_travel, na.rm = TRUE),
    median_time_travel = median(time_travel, na.rm = TRUE)
  )

endline.data %>%
  group_by(dist.pot.group) %>%
  summarise(
    frac_pay_0 = mean(travel_pay == 0, na.rm = TRUE),
    mean_pay = mean(travel_pay, na.rm = TRUE)
  )


endline.data = endline.data %>%
  mutate(
    travel_clean = case_when(
      travel == "1" ~ "walk",
      travel == "2" ~ "motorbike",
      travel == "3" ~ "car/taxi",
      travel == "4" ~ "bus",
      travel == "5" ~ "free ride"
    )
  )

  endline.data %>%
    write_csv("temp-data/endline-data.csv")

endline.data %>%
  count(travel_clean) %>%
  mutate(
    pct = 100*n/sum(n)
  )

endline.data %>%
  summarise(
    mean_time_travel = mean(time_travel, na.rm = TRUE),
    median_time_travel = median(time_travel, na.rm = TRUE)
  )


endline.data %>%
  summarise(
    frac_pay_0 = mean(travel_pay == 0, na.rm = TRUE),
    mean_pay = mean(travel_pay, na.rm = TRUE)
  )

endline.data %>%
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

analysis_data = analysis_data %>%
  left_join(
    pretreat_data %>%
      select(KEY.individ, 
      all_of(pretreat_vars)
      ),
      by = "KEY.individ"
  )

#### Functions for RF bootstrap ------------------------------------------------
# Estimate and Predict to generate ATEs

# For a given set of IDs, create bs data - n.b. this allows cluster to appear 
# multiple times
create_bs_data = function(split_data, ids) {
    bs_df = do.call(rbind, split_data[ids])
    return(bs_df)
}


# Create bootstrapped predictions
create_bs_preds = function(pred_df, ...) {
    # check if combined already present for PEAs
    combined_present = "combined" %in% pred_df$assigned_dist_group
    preds = pred_df %>%
        group_by(
            assigned_treatment,
            assigned_dist_group,
            ...
        ) %>%
        summarise(
            mean_pred = mean(pred, na.rm = TRUE),
            .groups = "drop"
        )
    if (!combined_present) {
      preds = preds %>%
        bind_rows(
            pred_df %>%
                group_by(
                    assigned_treatment,
                    ...
                ) %>%
                summarise(
                    mean_pred = mean(pred, na.rm = TRUE),
                    .groups = "drop"
                ) %>%
                mutate(
                    assigned_dist_group = "combined"
                )
        )
    }
        

    signal_pred = pred_df %>%
        mutate(
          signal = case_when(
            assigned_treatment %in% c("bracelet", "ink") ~ "signal",
            assigned_treatment %in% c("calendar", "control") ~ "no signal"
          )
        ) %>%
        group_by(signal, assigned_dist_group, ...) %>%
        summarise(
            mean_pred = mean(pred, na.rm = TRUE),
            .groups = "drop"
        ) 
            
    if (!combined_present) {
      signal_pred = signal_pred %>%
            bind_rows(
            pred_df %>%
                mutate(
                  signal = case_when(
                    assigned_treatment %in% c("bracelet", "ink") ~ "signal",
                    assigned_treatment %in% c("calendar", "control") ~ "no signal"
                  )
                ) %>%
                group_by(
                    signal, 
                    ...
                ) %>%
                summarise(
                    mean_pred = mean(pred, na.rm = TRUE)
                ) %>%
                mutate(
                    assigned_dist_group = "combined"
                )
            )
    } 
          
    preds = bind_rows(preds, signal_pred)
       return(preds)
}



# Function to generate samples from a Dirichlet distribution
generate_dirichlet <- function(alpha, n) {
  # Generate matrix of gamma samples
  gamma_samples <- matrix(rgamma(n * length(alpha), shape = alpha, scale = 1), nrow = n)
  
  # Normalize rows to sum to 1
  dirichlet_samples <- gamma_samples / rowSums(gamma_samples)
  
  return(dirichlet_samples)
}

pred_bs_f = function(f, data, weights, realised_fit = FALSE) {
    if (realised_fit == TRUE) {
        data$wt = 1
    } else {
        data$wt = weights[data$cluster_id]
    }
    fit = f(data, weights = ~wt)

    pred_data = bind_rows(
      data %>%
        mutate(assigned_treatment = "bracelet"),
      data %>%
        mutate(assigned_treatment = "calendar"),
      data %>%
        mutate(assigned_treatment = "ink"),
      data %>%
        mutate(assigned_treatment = "control")
    )

    pred_data$pred = predict(fit, newdata = pred_data)
    data = pred_data %>%
        select(
            assigned_dist_group,
            assigned_treatment,
            standard_cluster.dist.to.pot,
            pred,
            any_of("sms_treatment")
        )


    return(data)
}

pred_bs_f_at_x = function(f,  data, weights, realised_fit = FALSE) {
    if (realised_fit == TRUE) {
        data$wt = 1
    } else {
        data$wt = weights[data$cluster_id]
    }
    fit = f(data, weights = ~wt)
    data = data %>%
        group_by(assigned_dist_group) %>%
        mutate(
          average_standard_cluster.dist.to.pot = mean(standard_cluster.dist.to.pot),
          average_dist.to.pot = mean(dist.to.pot)
        ) %>%
        ungroup() 


    collapsed_data = data %>%  
      mutate(
        standard_cluster.dist.to.pot = average_standard_cluster.dist.to.pot,
        dist.to.pot = average_dist.to.pot,
        county = "Kakamega"
      ) %>%
      select(
        standard_cluster.dist.to.pot,
        dist.to.pot,
        county,
        assigned_dist_group,
        assigned_treatment,
        signal,
        any_of("sms_treatment")
      ) %>%
      distinct() %>%
      bind_rows(
        data %>%
          group_by(assigned_treatment) %>%
          summarise(
              standard_cluster.dist.to.pot = mean(standard_cluster.dist.to.pot),
              dist.to.pot = mean(dist.to.pot),
              county = "Kakamega",
              signal = unique(signal)
          ) %>%
          mutate(
            assigned_dist_group = "combined"
          )
      )

    pea_pred = collapsed_data %>%  
      predict(
        fit, newdata = .
      )

    collapsed_data$pred = pea_pred


    pea_signal_pred = collapsed_data %>%
      predict(
        fit, newdata = .
      )

    collapsed_data$signal_pred = pea_signal_pred
    collapsed_data = collapsed_data %>%
        select(
            assigned_dist_group,
            assigned_treatment,
            signal,
            standard_cluster.dist.to.pot,
            pred,
            signal_pred,
            any_of("sms_treatment")
        )

    return(collapsed_data)

}

bayes_bs_f_at_x = function(seed, f, data, ...) {
    set.seed(seed)
    n_clusters = length(unique(data$cluster.id))
    alpha = rep(1, n_clusters)
    weights = generate_dirichlet(alpha, 1)
    bs_fit = pred_bs_f_at_x(f,  data, weights = weights) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 

bayes_bs_f = function(seed, f, data, ...) {
    set.seed(seed)
    n_clusters = length(unique(data$cluster.id))
    alpha = rep(1, n_clusters)
    weights = generate_dirichlet(alpha, 1)
    bs_fit = pred_bs_f(f, data, weights = weights) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 

actual_bayesian_bs_fit_at_x = function(seed, f,  data, ...) {
    bs_fit = pred_bs_f_at_x(f,  data, 1, realised_fit = TRUE) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 
actual_bayesian_bs_fit = function(seed, f,  data, ...) {
    bs_fit = pred_bs_f(f,  data, 1, realised_fit = TRUE) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    data$wt = 1
    fit = f(data, weights = ~wt)
    return(list(bs_fit = bs_fit, fit = fit))
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
          conf.low = quantile(estimate, (1 - ci_width)/2, na.rm = TRUE),
          conf.high = quantile(estimate, 1 - (1 - ci_width)/2, na.rm = TRUE)
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
create_regression_output = function(data, f,  B_draws = 500, 
                                    stat = params$stat,
                                    caption = "Average Treatment Effects: Reduced Form",
                                    dependent_var = "Dependent variable: Take-up",
                                    type = "APE",
                                    stars = TRUE,
                                    drop_H0s = FALSE
                                    ) {

  if (type == "APE") {
    bs_f = bayes_bs_f
    actual_f = actual_bayesian_bs_fit
  } else {
    bs_f = bayes_bs_f_at_x
    actual_f = actual_bayesian_bs_fit_at_x
  }
  bs_draws = map_dfr(
    1:B_draws,
    ~bs_f(
      seed = .x,
      f = f,
      data = data
    ),
    .progress = TRUE
    )

  clean_te_draws_df = bs_draws %>%
    clean_te_draws()

  clean_signal_draws_df = bs_draws %>%
    clean_signal_draws()

  realised_fit_output = actual_f(
    seed = "realised fit",
    f = f,
    data = data
  )
  realised_fit = realised_fit_output$bs_fit
  n_obs = nobs(realised_fit_output$fit)



  signal_fit = realised_fit %>%
    clean_signal_draws() %>%
    rename(realised_pred = estimate) %>%
    select(realised_pred, assigned_dist_group, assigned_treatment)


  te_fit = realised_fit %>%
    clean_te_draws() %>%
    rename(realised_pred = estimate) %>%
    select(realised_pred, assigned_dist_group, assigned_treatment)

  signal_summ = add_summ_stats(clean_signal_draws_df, signal_fit)
  te_summ = add_summ_stats(clean_te_draws_df, te_fit)

  pval_only_terms = c("bracelet - calendar", "signal")

  n_obs_df = tibble(
    assigned_treatment = factor("Observations"),
    pval = as.character(prettyNum(n_obs, big.mark = ",")),
    assigned_dist_group = unique(te_summ$assigned_dist_group),
    show_pval_only = TRUE,
    n_obs_line = TRUE
  )

  overall_summ = bind_rows(
    signal_summ,
    te_summ
  ) %>%
    mutate(
      show_pval_only = assigned_treatment %in% pval_only_terms
    ) %>%
    filter(assigned_treatment != "no signal")  %>%
    mutate(n_obs_line = FALSE) %>%
    bind_rows(n_obs_df) 

  if (drop_H0s) {
    overall_summ = overall_summ %>%
      filter(show_pval_only == FALSE | n_obs_line == TRUE)
  }

  default_tbl = overall_summ %>%
    prep_tbl(stat = stat, stars = stars) %>%
    nice_kbl_table(
      cap = caption,
      outcome_var = dependent_var
      )


  different_order_tbl = overall_summ %>%
    prep_tbl(stat = stat, stars = stars) %>%
    filter(assigned_treatment != "Observations") %>%
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
          "Ink", "Bracelet", "Calendar"
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
prep_tbl = function(tes, stat = "ci", stars = FALSE) {

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
        "Observations",
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
            mutate(val = paste0("{(", round(std_error, 3), ")}"))
    } else {
        tes = tes %>%
            mutate(
                val = paste0("{(", round_pval(pval, 3), ")}")
            )
    }




    tbl =  tes %>%
        select(assigned_treatment, assigned_dist_group, estimate, conf.low, conf.high, val, pval, oneside_pval, show_pval_only, n_obs_line)  %>%
        mutate(
          show_stars = 
            ((assigned_treatment %in% c("bracelet", "calendar", "ink")) | assigned_dist_group == "far - close") & stars == TRUE,
          stars = case_when(
            pval < 0.001 ~ "***",
            pval < 0.05 ~ "**",
            pval < 0.1 ~ "*",
            TRUE ~ ""
          )
        ) %>%
        mutate(across(where(is.numeric), ~round(.x, 3))) %>%
        mutate(
            estim_std = case_when(
              show_pval_only == TRUE ~ linebreak(paste0(pval), align = "c"),
              show_stars == TRUE ~ linebreak(paste0(estimate, stars, "\n", str_glue("{val}")), align = "c"),
              TRUE ~ linebreak(paste0(estimate,"\n", str_glue("{val}")), align = "c")
            )
        ) %>%
        bind_rows(
            filter(., show_pval_only == TRUE & n_obs_line == FALSE) %>%
              mutate(
                estim_std = linebreak(paste0(pval), align = "c")
              ) %>%
              mutate(estim_std = as.character(estim_std)) %>%
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
        )  %>%
        filter(
          !(assigned_treatment %in% c("Bracelet - No Signal", "$H0$: Bracelet > Calendar, $p$-value", "$H0$: Any Signal > No Signal, $p$-value"))
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
  row_spec(c(3), hline_after = TRUE) 
}

custom_save_latex_table = function(table, table_name, table_output_path = params$table_output_path){
  table_conn = file(
    file.path(
      table_output_path, paste0(table_name, ".tex")
    )
  )

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

wrapper_function = function(data, regression_spec, tidy_summ_path, table_name, table_options = list()) {
  default_table_options = list(
    caption = "Average Treatment Effects: Reduced Form",
    dependent_var = "Dependent variable: Take-up",
    type = "APE",
    stars = TRUE,
    drop_H0s = FALSE
  )
  table_options = modifyList(default_table_options, table_options)
  output = create_regression_output(
    data = data,
    f = regression_spec,
    caption = table_options$caption,
    dependent_var = table_options$dependent_var,
    type = table_options$type,
    stars = table_options$stars,
    drop_H0s = table_options$drop_H0s

  )
  output$tidy_summary %>%
    write_csv(tidy_summ_path)
  output$default_tbl %>%
    custom_save_latex_table(
      table_name = table_name
    )
  output$different_order_tbl %>%
    custom_save_latex_table(
      table_name = paste0(table_name, "_weird_order")
    )
    return(output)
}

#### Frequentist Estimates ####
# load expected distances
cell_expected_dist_df = read_csv(here::here("data", "cell_expected_dist_df.csv")) %>%
  rename(cell_expected_dist = dist, assigned_treatment = random_treat, assigned_dist_group = random_dist_group) 
cluster_expected_dist_df = read_csv(here::here("data", "cluster_expected_dist.csv")) %>%
  rename(clust_expected_dist = dist) %>%
  mutate(cluster.id = factor(cluster.id))


analysis_data = analysis_data %>%
    mutate(
        county = factor(county),
        cluster.id = factor(cluster.id),
        assigned_treatment = assigned.treatment,
        assigned_dist_group = dist.pot.group,
        signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
        signal = factor(signal, levels = c("no signal", "signal"))
    ) %>%
    left_join(
        cluster_expected_dist_df,
        by = "cluster.id"
    ) %>%
    left_join(
        cell_expected_dist_df,
        by = c("assigned_treatment", "assigned_dist_group")
    ) %>%
    mutate(
      standard_clust_expected_dist = clust_expected_dist/sd_of_dist,
      standard_cell_expected_dist = cell_expected_dist/sd_of_dist
    )
analysis_data %>%
  write_csv(
    "temp-data/analysis-data.csv"
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

## Covariates Added
cov_analysis_data = read_csv("temp-data/analysis-cluster-covariate-data.csv") %>%
  mutate(assigned_dist_group = dist.pot.group) %>%
  mutate(
    cluster.id = cluster_id
  ) %>%
  mutate(assigned.treatment = factor(assigned.treatment, levels = c("control", "ink", "calendar", "bracelet"))) %>%
  mutate(assigned_treatment = assigned.treatment)  %>%
  mutate(wt = 1) %>%
  mutate(
      county = factor(county),
      cluster.id = factor(cluster.id),
      assigned_treatment = assigned.treatment,
      assigned_dist_group = dist.pot.group,
      signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
      signal = factor(signal, levels = c("no signal", "signal"))
  ) %>%
  left_join(
      cluster_expected_dist_df %>%
        mutate(cluster.id = as.numeric(cluster.id)),
      by = c("cluster_id" = "cluster.id")
  ) %>%
  left_join(
      cell_expected_dist_df,
      by = c("assigned_treatment", "assigned_dist_group")
  ) %>%
  mutate(
    standard_clust_expected_dist = clust_expected_dist/sd_of_dist,
    standard_cell_expected_dist = cell_expected_dist/sd_of_dist
  ) %>%
  mutate(
    mu_d = standard_clust_expected_dist
  )

write_csv(cov_analysis_data, "temp-data/analysis-cluster-recentered-covariate-data.csv")

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

#### Takeup Continuous Distance + LASSO Covs + Cluster Expected Distance
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


#### Beliefs -------------------------------------------------------------------
belief_ana_df = analysis_data %>%
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
    filter(obs_know_person > 0)


all_data = analysis.data %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


disagg_belief_all_df = all_data %>%
  mutate(
    female = gender == "female"
  ) %>%
  left_join(
    cov_analysis_data %>%
      select(cluster.id,  mu_d, standard_clust_expected_dist) %>%
      mutate(cluster.id = as.numeric(cluster.id)) %>%
      unique(),
    by = "cluster.id"
  ) %>% 
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
      county,
      dewormed,
      all_of(l_cov_vars),
      standard_clust_expected_dist
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
           county,
           dewormed,
           all_of(l_cov_vars),
           standard_clust_expected_dist
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
    


disagg_base_belief_data = cov_analysis_data %>%
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
      county,
      dewormed,
      all_of(l_cov_vars),
      standard_clust_expected_dist
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
           county,
           dewormed,
           all_of(l_cov_vars),
           standard_clust_expected_dist
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
  )  %>%
  left_join(
      cluster_expected_dist_df %>%
        mutate(cluster.id = as.numeric(cluster.id)),
      by = c("cluster_id" = "cluster.id")
  ) %>%
    mutate(
      standard_clust_expected_dist = clust_expected_dist/sd_of_dist,
      mu_d = standard_clust_expected_dist
    )


know_all_df = disagg_belief_all_df %>%
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
  )  %>%
  left_join(
      cluster_expected_dist_df %>%
        mutate(cluster.id = as.numeric(cluster.id)),
      by = c("cluster_id" = "cluster.id")
  ) %>%
    mutate(
      standard_clust_expected_dist = clust_expected_dist/sd_of_dist,
      mu_d = standard_clust_expected_dist
    )


know_1_df = know_df  %>%
  filter(belief_type == "1ord") 
know_2_df = know_df  %>%
  filter(belief_type == "2ord")


know_1_all_df = know_all_df  %>%
  filter(belief_type == "1ord") %>%
  mutate(cluster_id = as.numeric(cluster.id)) 

know_1_all_df %>%
  select(cluster_id)

discrete_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] +  mu_d | county,
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
discrete_fob_output = wrapper_function(
  data = know_df %>%
    filter(belief_type == "1ord"),
  regression_spec = discrete_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_discrete_fob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-fob-tidy-tes.csv"
)

#### Checking using full sample (SMS + non-monitored) doesn't change results

discrete_fob_full_output = wrapper_function(
  data = know_1_all_df,
  regression_spec = discrete_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_discrete_fob_fullsample_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-fob-fullsample-tidy-tes.csv"
)

#### Checking weird 150 individuals don't cause issues

not_in_monitored = readRDS("temp-data/not_in_monitored.rds")

clean_endline_extra_df =  endline.know.table.data %>% 
    filter(fct_match(know.table.type, "table.A")) %>%
    filter(KEY.individ %in% not_in_monitored)  %>%
    group_by(KEY.individ) %>%
    summarise(
      obs_know_person = sum(num.recognized),
      obs_know_person_prop = mean(num.recognized),
      knows_other_dewormed = sum(fct_match(dewormed, c("yes", "no")), na.rm = TRUE),
      knows_other_dewormed_yes = sum(fct_match(dewormed, "yes"), na.rm = TRUE),
      knows_other_dewormed_no = sum(fct_match(dewormed, "no"), na.rm = TRUE),
      thinks_other_knows = sum(fct_match(second.order, c("yes", "no")), na.rm = TRUE),
      thinks_other_knows_yes = sum(fct_match(second.order, "yes"), na.rm = TRUE),
      thinks_other_knows_no = sum(fct_match(second.order, "no"), na.rm = TRUE),
      assigned.treatment = first(assigned.treatment),
      dist.pot.group = first(dist.pot.group),
      cluster.id = first(cluster.id)
    ) %>%
    left_join(
      cluster_expected_dist_df %>%
        mutate(cluster.id = as.numeric(cluster.id)),
      by = c("cluster.id" = "cluster.id")
    ) %>%
    # convert to same format as know_df
    mutate(
      standard_clust_expected_dist = clust_expected_dist/sd_of_dist,
      mu_d = standard_clust_expected_dist,
    ) %>%
    mutate(
      assigned_treatment = assigned.treatment,
      assigned_dist_group = dist.pot.group,
      signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
      signal = factor(signal, levels = c("no signal", "signal"))
    ) %>%
    mutate(
      prop_knows = obs_know_person_prop,
      cluster.id = factor(cluster.id)
    ) %>%
    left_join(
      analysis_data %>%
        select(cluster.id, county, standard_cluster.dist.to.pot) %>%
        unique(), 
        by = "cluster.id"
    )



extra_know_df = bind_rows(
  know_df %>% 
    filter(belief_type == "1ord"),
  clean_endline_extra_df %>%
    mutate(cluster.id = factor(cluster.id))
)


extra_fit = extra_know_df %>%
  feols(
    prop_knows ~ assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control")  ,
    cluster = ~cluster.id
  )
single_fit = know_df %>%
  filter(belief_type == "1ord") %>%
  feols(
    prop_knows ~ assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") ,
    cluster = ~cluster.id
  )


tidy_comp_extra_df = bind_rows(
  extra_fit %>%
    tidy(conf.int = TRUE) %>%
    mutate(type = "extra 150 people"),
    single_fit %>%
    tidy(conf.int = TRUE) %>% 
    mutate(type = "current analysis")
)

tidy_comp_extra_df %>%
  ggplot(aes(
    x = estimate,
    xmin = conf.low,
    xmax = conf.high,
    y = term,
    colour = type
  )) +
  geom_pointrange(
    position = position_dodge(width = 0.5),
    size = 1
  )  +
  theme_bw()
ggsave(
  "temp-data/extra_fit_comp.pdf",
  width = 10,
  height = 5
)



#### FOB Continuous Distance + LASSO Covs + Cluster Expected Distance
cts_fob_output = wrapper_function(
  data = know_df %>%
    filter(belief_type == "1ord"),
  regression_spec = cts_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_cts_fob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-cts-fob-tidy-tes.csv"
)

#### FOB HH Distance + LASSO Covs + Cluster Expected Distance
hh_fob_output = wrapper_function(
  data = know_df %>%
    filter(belief_type == "1ord"),
  regression_spec = hh_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability"
  ),
  table_name = "rf_hh_fob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-hh-fob-tidy-tes.csv"
)

#### SOB Discrete Distance + LASSO Covs + Cluster Expected Distance
discrete_sob_output = wrapper_function(
  data = know_df %>%
    filter(belief_type == "2ord"),
  regression_spec = discrete_f_know,
  table_options = list(
    dependent_var = "Dependent variable: Observability Beliefs"
  ),
  table_name = "rf_discrete_sob_spec_tbl",
  tidy_summ_path = "temp-data/tidy-rf-tes/reducedform-discrete-sob-tidy-tes.csv"
)



#### Levels --------------------------------------------------------------------

#### Takeup LEVELS Continuous Distance + LASSO Covs + Expected Distance
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



#### FOB Levels Discrete Distance + LASSO Covs + Expected Distance
fob_bs_draws = map_dfr(
  1:500,
  ~bayes_bs_f(
    seed = .x,
    f = discrete_f_know,
    data = know_df %>%
      filter(belief_type == "1ord")
  ),
  .progress = TRUE
  )

fob_levels_point = actual_bayesian_bs_fit(
  seed = "realised fit",
  f = discrete_f_know,
  data = know_df %>%
    filter(belief_type == "1ord")
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



#### Alternative Regressions ---------------------------------------------------

####  Endline Predicted Deworming Takeup
endline_data = endline.data %>%
  mutate(
    assigned_treatment = as_factor(assigned.treatment), 
    assigned_dist_group = as_factor(dist.pot.group),
    cluster_id = as_factor(cluster.id),
    # this isn't actually used
    standard_cluster.dist.to.pot = dist.to.pot/sd_of_dist,
    dworm_frac = dworm_rate / 10,
    # different naming convention here
    have_ink = ink_visible
  )  %>%
  left_join(
    cov_analysis_data %>%
      select(KEY.individ, mu_d, all_of(l_cov_vars)),
      by = "KEY.individ"
  )



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

long_incentive_check_df = endline_data %>%
  select(all_of(c(got_vars, have_vars, seen_vars, mean_vars)), assigned_treatment, cluster_id, county)  %>%
  pivot_longer(
    cols = all_of(c(got_vars, have_vars, seen_vars, mean_vars))
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
    got, # "Did you receive X when you went for deworming?"
    have, # "Do you still have X?"
    seen, # "have you seen people wearing these "X"?"/ "have you seen people these calendars before"
    meandeworm # "What does it mean if a person has a bracelet?" -> coded into deworm mentions
  )



n_incentive_df %>%
  pivot_wider(
    names_from = variable_type,
    values_from = n
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
        "Received incentive when treated", 
        "Have incentive currently", 
        "Seen incentive", 
        "Link incentive to deworming"),
      escape = FALSE, 
      booktabs = TRUE,
      align = "lcccc", 
      caption = "Endline Incentive Checks"
  ) %>% 
  row_spec(c(2), hline_after = TRUE) 


incentive_check_tbl %>%
  custom_save_latex_table("incentive-check-tbl")

#### Preference for Gift Fit Not Dewormed ---------------------------------------
pref_gift_fit_not_dewormed = analysis_data %>%
    filter(!is.na(gift_choice), monitored, monitor.consent, !hh.baseline.sample.pool, !is.na(sms.treatment)) %>% 
    group_by(assigned.treatment, dist.pot.group, dewormed) %>% 
    mutate(arm.size = n()) %>% 
    group_by(gift_choice, add = TRUE) %>%
    # filter(assigned.treatment %in% c("control",  "calendar", "bracelet")) %>%
    filter(
      dewormed == FALSE
    )  %>%
    ungroup() %>%
    select(KEY.individ, cluster.id, gift_choice, assigned.treatment, dist.pot.group, county, standard_cluster.dist.to.pot) %>%
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


wrapper_function(
  data = pref_gift_fit_not_dewormed,
  regression_spec = pref_gift_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/preference-for-bracelet-tidy-tes.csv",
  table_name = "preference_for_bracelet_spec_tbl",
  table_options = list(
    caption = "Average Treatment Effects: Reduced Form", 
    dependent_var = "Dependent variable: Prefer Bracelet", 
    stars = TRUE, 
    type = "APE",
    drop_H0s = TRUE
    )
)

#### Travel Time --------------------------------------------------------

endline.data %>%
  count(travel)

  endline.data %>%
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

travel_time_df = endline.data %>%
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

endline.data = endline.data %>%
  mutate(
    travel_clean = case_when(
      travel == "1" ~ "walk",
      travel == "2" ~ "motorbike",
      travel == "3" ~ "car/taxi",
      travel == "4" ~ "bus",
      travel == "5" ~ "free ride"
    )
  )


endline.data %>%
  count(travel_clean) %>%
  mutate(
    pct = 100*n/sum(n)
  )

endline.data %>%
  summarise(
    mean_time_travel = mean(time_travel, na.rm = TRUE),
    median_time_travel = median(time_travel, na.rm = TRUE)
  )


endline.data %>%
  summarise(
    frac_pay_0 = mean(travel_pay == 0, na.rm = TRUE),
    mean_pay = mean(travel_pay, na.rm = TRUE)
  )

endline.data %>%
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

ggsave("temp-data/p-sms-tes.pdf", width = 8, height = 6)

#### Heterogeneity -------------------------------------------------------------
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
  )

analysis_data = analysis_data %>%
  mutate(
    age_gt_40 = age.census > 40
  ) %>%
  mutate(treatment = factor(assigned_treatment, levels = c("control", "bracelet", "calendar", "ink"))) 

clean_perception_data = baseline.data %>% 
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


  


age_het_fit = analysis_data %>%
  feols(
    dewormed ~  
      treatment + 
      standard_cluster.dist.to.pot + 
      age_gt_40 +
      i(treatment, age_gt_40, "control")  
      | county,
      cluster = ~cluster.id
  ) 


judge_het_fit = analysis_data %>%
  left_join(
    overall_judgement_score_df %>% 
      select(cluster.id, judge_score_gt_mean),
      by = "cluster.id"
  ) %>%
  feols(
    dewormed ~  
      treatment + 
      standard_cluster.dist.to.pot + 
      judge_score_gt_mean +
      i(treatment, judge_score_gt_mean, "control")  
      | county,
      cluster = ~cluster.id
  ) 

phone_het_fit = analysis_data %>%
  feols(
    dewormed ~  
      treatment + 
      standard_cluster.dist.to.pot + 
      have_phone_lgl +
      i(treatment, have_phone_lgl, "control")  
      | county,
      cluster = ~cluster.id
  ) 
  


prevdeworm_het_fit = analysis_data %>%
  feols(
    dewormed ~  
      treatment + 
      standard_cluster.dist.to.pot + 
      frac_prev_dewormed_gt_mean +
      i(treatment, frac_prev_dewormed_gt_mean, "control")  
      | county,
      cluster = ~cluster.id
  )
  
externality_het_fit = analysis_data %>%
  feols(
    dewormed ~  
      treatment + 
      standard_cluster.dist.to.pot + 
      frac_externality_gt_mean +
      i(treatment, frac_externality_gt_mean, "control")  
      | county,
      cluster = ~cluster.id
  )


gender_het_fit = analysis_data %>%
  mutate(female = gender == "female") %>%
  feols(
    dewormed ~  
      treatment + 
      standard_cluster.dist.to.pot + 
      female +
      i(treatment, female, "control")  
      | county,
      cluster = ~cluster.id
  )

tex_postprocessing = function(tex) {
    tex %>%
        str_remove("\\\\begin\\{table\\}\\[htbp\\]") %>%
        str_remove("\\\\end\\{table\\}") %>%
        str_replace(
          .,
          "Covariate",
          "\\\\midrule Covariate"
        )
}

  etable(
    list(
      "Age $>40$" = age_het_fit,
      "Female" = gender_het_fit,
      "Phone Owner" = phone_het_fit,
      "Judgement" = judge_het_fit,
      "Previously Dewormed" = prevdeworm_het_fit,
      "Externality Knowledge" = externality_het_fit
    ),
    drop = "dist",
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
    postprocess.tex = tex_postprocessing,
    replace = TRUE,
    style.df = style.df(
      depvar.title = "", 
      fixef.title = "")
  )

feols(
  data = analysis_data,
  dewormed ~ treatment*standard_cluster.dist.to.pot  + age_gt_40 | county,
  cluster = ~cluster.id
)


age_het_preds =  bind_rows(
  age_het_fit %>%
    avg_comparisons(
      newdata = datagrid(
        age_gt_40 =  FALSE
      ),
      variables = list(
        assigned_treatment = "reference"
      )
    ) %>%
    tidy(conf.int = TRUE),
  age_het_fit %>%
    avg_predictions(
      newdata = datagrid(
        assigned_treatment = "control"
      ),
      variables = list(age_gt_40 = c(TRUE, FALSE))
    ) %>%
    tidy(conf.int = TRUE)  %>%
    mutate(
      control_mean = paste0("control_", age_gt_40)
    ),
  age_het_fit %>%
    avg_comparisons(
      variables = list(
        assigned_treatment = "reference",
        age_gt_40  = c(TRUE, FALSE)
        ),
      cross = TRUE
    )  %>%
    tidy(
      conf.int = TRUE,
    )
) %>%
    mutate(
      het_variable = "age_gt_40"
    )




judge_het_fit = analysis_data %>%
  left_join(
    overall_judgement_score_df %>% 
      select(cluster.id, judge_score_gt_mean),
      by = "cluster.id"
  ) %>%
  feglm(
    dewormed ~ 0 + 
      assigned_treatment + 
      standard_cluster.dist.to.pot + 
      i(assigned_treatment, standard_cluster.dist.to.pot, "control") +
      judge_score_gt_mean +
      i(assigned_treatment, judge_score_gt_mean, "control")  
      | county,
      family = "probit",
      cluster = ~cluster.id
  ) 
  

judge_het_preds =  bind_rows(
  judge_het_fit %>%
    avg_comparisons(
      newdata = datagrid(
        judge_score_gt_mean =  FALSE
      ),
      variables = list(
        assigned_treatment = "reference"
      )
    ) %>%
    tidy(conf.int = TRUE),
  judge_het_fit %>%
    avg_predictions(
      newdata = datagrid(
        assigned_treatment = "control"
      ),
      variables = list(judge_score_gt_mean = c(TRUE, FALSE))
    ) %>%
    tidy(conf.int = TRUE)  %>%
    mutate(
      control_mean = paste0("control_", judge_score_gt_mean)
    ),
  judge_het_fit %>%
    avg_comparisons(
      variables = list(
        assigned_treatment = "reference",
        judge_score_gt_mean  = c(TRUE, FALSE)
        ),
      cross = TRUE
    )  %>%
    tidy(
      conf.int = TRUE,
    )
) %>%
    mutate(
      het_variable = "judge_score_gt_mean"
    )




phone_het_fit = analysis_data %>%
  feglm(
    dewormed ~ 0 + 
      assigned_treatment + 
      standard_cluster.dist.to.pot + 
      i(assigned_treatment, standard_cluster.dist.to.pot, "control") +
      have_phone_lgl +
      i(assigned_treatment, have_phone_lgl, "control")  
      | county,
      family = "probit",
      cluster = ~cluster.id
  ) 
  
  
  
phone_het_preds =  bind_rows(
  phone_het_fit %>%
    avg_comparisons(
      newdata = datagrid(
        have_phone_lgl =  FALSE
      ),
      variables = list(
        assigned_treatment = "reference"
      )
    ) %>%
    tidy(conf.int = TRUE),
  phone_het_fit %>%
    avg_predictions(
      newdata = datagrid(
        assigned_treatment = "control"
      ),
      variables = list(have_phone_lgl = c(TRUE, FALSE))
    ) %>%
    tidy(conf.int = TRUE)  %>%
    mutate(
      control_mean = paste0("control_", have_phone_lgl)
    ),
  phone_het_fit %>%
    avg_comparisons(
      variables = list(
        assigned_treatment = "reference",
        have_phone_lgl  = c(TRUE, FALSE)
        ),
      cross = TRUE
    )  %>%
    tidy(
      conf.int = TRUE,
    )
) %>%
    mutate(
      het_variable = "have_phone_lgl"
    )


prevdeworm_het_fit = analysis_data %>%
  feglm(
    dewormed ~ 0 + 
      assigned_treatment + 
      standard_cluster.dist.to.pot + 
      i(assigned_treatment, standard_cluster.dist.to.pot, "control") +
      frac_prev_dewormed_gt_mean +
      i(assigned_treatment, frac_prev_dewormed_gt_mean, "control")  
      | county,
      family = "probit",
      cluster = ~cluster.id
  )
  
prevdeworm_het_preds =  bind_rows(
   prevdeworm_het_fit %>%
    avg_comparisons(
      newdata = datagrid(
        frac_prev_dewormed_gt_mean =  FALSE
      ),
      variables = list(
        assigned_treatment = "reference"
      )
    ) %>%
    tidy(conf.int = TRUE),
  prevdeworm_het_fit %>%
    avg_predictions(
      newdata = datagrid(
        assigned_treatment = "control"
      ),
      variables = list(frac_prev_dewormed_gt_mean = c(TRUE, FALSE))
    ) %>%
    tidy(conf.int = TRUE)  %>%
    mutate(
      control_mean = paste0("control_", frac_prev_dewormed_gt_mean)
    ),
  prevdeworm_het_fit %>%
    avg_comparisons(
      variables = list(
        assigned_treatment = "reference",
        frac_prev_dewormed_gt_mean  = c(TRUE, FALSE)
        ),
      cross = TRUE
    )  %>%
    tidy(
      conf.int = TRUE,
    )
) %>%
    mutate(
      het_variable = "frac_prev_dewormed_gt_mean"
    )
  
  

externality_het_fit = analysis_data %>%
  feglm(
    dewormed ~ 0 + 
      assigned_treatment + 
      standard_cluster.dist.to.pot + 
      i(assigned_treatment, standard_cluster.dist.to.pot, "control") +
      frac_externality_gt_mean +
      i(assigned_treatment, frac_externality_gt_mean, "control")  
      | county,
      family = "probit",
      cluster = ~cluster.id
  )
  
  
  
  
  
externality_het_preds =  bind_rows(
  externality_het_fit %>%
    avg_comparisons(
      newdata = datagrid(
        frac_externality_gt_mean =  FALSE
      ),
      variables = list(
        assigned_treatment = "reference"
      )
    ) %>%
    tidy(conf.int = TRUE),
  externality_het_fit %>%
    avg_predictions(
      newdata = datagrid(
        assigned_treatment = "control"
      ),
      variables = list(frac_externality_gt_mean = c(TRUE, FALSE))
    ) %>%
    tidy(conf.int = TRUE)  %>%
    mutate(
      control_mean = paste0("control_", frac_externality_gt_mean)
    ),
  externality_het_fit %>%
    avg_comparisons(
      variables = list(
        assigned_treatment = "reference",
        frac_externality_gt_mean  = c(TRUE, FALSE)
        ),
      cross = TRUE
    )  %>%
    tidy(
      conf.int = TRUE,
    )
) %>%
    mutate(
      het_variable = "frac_externality_gt_mean"
    )
  

gender_het_fit = analysis_data %>%
  feglm(
    dewormed ~ 0 + 
      assigned_treatment + 
      standard_cluster.dist.to.pot + 
      i(assigned_treatment, standard_cluster.dist.to.pot, "control") +
      gender +
      i(assigned_treatment, gender, "control")  
      | county,
      family = "probit",
      cluster = ~cluster.id
  )
    
gender_het_preds =  bind_rows(
  gender_het_fit %>%
    avg_comparisons(
      newdata = datagrid(
        gender = "male"
      ),
      variables = list(
        assigned_treatment = "reference"
      )
    ) %>%
    tidy(conf.int = TRUE),
  gender_het_fit %>%
    avg_predictions(
      newdata = datagrid(
        assigned_treatment = "control"
      ),
      variables = list(gender = c("male", "female"))
    ) %>%
    tidy(conf.int = TRUE)  %>%
    mutate(
      control_mean = paste0("control_", gender)
    ),
  gender_het_fit %>%
    avg_comparisons(
      variables = list(
        assigned_treatment = "reference",
        gender  = c("male", "female")
        ),
      cross = TRUE
    )  %>%
    tidy(
      conf.int = TRUE,
    )
) %>%
    mutate(
      het_variable = "gender"
    )





het_preds = bind_rows(
  age_het_preds,
  gender_het_preds,
  phone_het_preds,
  judge_het_preds,
  prevdeworm_het_preds,
  externality_het_preds
) %>%
  select(het_variable, control_mean, contrast, assigned_treatment, contrast_assigned_treatment, estimate, std.error, p.value )


het_tbl = het_preds %>%
  mutate(
    te_diff = !is.na(contrast_assigned_treatment),
    te_diff = case_when(
      control_mean == "control_female" ~ TRUE,
      control_mean == "control_male" ~ FALSE,
      str_detect(control_mean, "TRUE") ~ TRUE,
      str_detect(control_mean, "FALSE") ~ FALSE,
      TRUE ~ te_diff
    ),
    term = case_when(
      !is.na(control_mean) ~ "Control mean",
      !is.na(contrast) & contrast != "control" ~ str_extract(contrast, "(?<=mean\\()\\w+") %>% str_to_title(),
      !is.na(contrast_assigned_treatment)  ~ str_extract(contrast_assigned_treatment, "(?<=mean\\()\\w+")  %>% str_to_title()
    )
  ) %>%
  select(het_variable, term, te_diff, estimate, std.error) %>%
  mutate(
    across(where(is.numeric), ~round(., 3)),
    val = paste0("{[", std.error, "]}"),
    estim_std = linebreak(paste0(estimate,"\n", str_glue("{val}")), align = "c") 
  )  %>%
  select(-estimate, -std.error, -val) %>%
  pivot_wider(
    names_from = c(het_variable, te_diff),
    values_from = estim_std
  ) %>%
  kbl(
    col.names = c(
      "Treatment",

      "$\\leq 40$",
      "$> 40$",

      "Male",
      "Female",

      "No",
      "Yes",

      "No",
      "Yes",

      "No",
      "Yes",

      "No",
      "Yes"
    ),
    align = "lcccccccccccc",
    booktabs = TRUE,
    format = "latex",
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("scale_down")
  ) %>% 
  add_header_above(
   c(
    " " = 1,
    "Age" = 2,
    "Gender" = 2,
    "Phone Owner" = 2,
    "Judgemental" = 2,
    "Previously Dewormed" = 2,
    "Understand Externalities" = 2
    )
  )  %>%
  add_header_above(
    c(
      " " = 7,
      "Community Level" = 6
    )
  ) %>%
  row_spec(c(3), hline_after = TRUE) 

het_tbl %>%
  custom_save_latex_table(
    table_name = "het-tes-tbl"
  )

het_tbl = het_fits %>%
  mutate(
    term = str_extract(
      contrast_assigned_treatment,
      "(?<=mean\\()\\w+"
    )
  ) %>%
  select(
    het_variable, term, estimate, std.error, p.value
  ) %>%
  mutate(
    across(where(is.numeric), ~round(., 3)),
    val = paste0("{[", std.error, "]}"),
    estim_std = linebreak(paste0(estimate,"\n", str_glue("{val}")), align = "c") 
  ) %>%
  select(het_variable, term, estim_std) %>%
  mutate(term = str_to_title(term))  %>%
  pivot_wider(
    names_from = het_variable,
    values_from = estim_std
  )  %>%
  kbl(
    col.names = c(
      "Heterogeneous Treatment Effects",
      "Female",
      "Phone Owner",
      "Community Judgemental",
      "Community Understand Externalities"
      ),
    booktabs = TRUE,
    escape = FALSE,
    align = "lcccc",
    format = "latex"
  ) %>%
  kable_styling(
    latex_options = c("scale_down")
  ) 



 het_tbl %>%
  custom_save_latex_table(
    table_name = "het-tes-tbl"
  )
