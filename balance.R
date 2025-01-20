#!/usr/bin/Rscript
script_options = docopt::docopt(
  stringr::str_glue("Usage:
  balance.R [options]
  
Options:
  --output-path=<path>  Where to save output files [default: {file.path('temp-data')}]
  --community-level
  --fit-ri
"),
  args = if (interactive()) "
    --output-path=temp-data \
    --fit-ri
    " else commandArgs(trailingOnly = TRUE)
) 


set.seed(12932)

library(tidyverse)
library(marginaleffects)
library(broom)
library(knitr)
library(kableExtra)
library(ggthemes)
library(fixest)
library(magrittr)
library(furrr)




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
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
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
# # creating a single treat x distance variable for balance testing
# cluster_treat_df = analysis_data %>%
#   mutate(
#       treat_dist = paste0(
#       "treat: ", 
#       assigned.treatment,
#       ", dist: ", dist.pot.group
#       ) %>% factor()
#   ) %>%
#   select(cluster.id, treat_dist, cluster.dist.to.pot) %>%
#   unique()


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


# Add church and add announce

# aggregate to household level (for messages Q)
# then to community
# first five (knowledge of deworming)
# last two (knowledge of availability)


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

# Adding school (PoT) data to analysis df
analysis_data = analysis_data %>%
    group_by(cluster.id) %>%
    mutate(row_id = 1:n()) %>%
    left_join(
      n_indiv_df %>% mutate(row_id = 1),
      by = c("cluster.id", "row_id")
    ) %>%
    ungroup() %>%
    select(-row_id) %>%
    mutate(have_phone_lgl = have_phone == "Yes")



analysis_data = analysis_data %>%
  mutate(
      treat_dist = paste0(
      "treat: ", 
      assigned.treatment,
      ", dist: ", dist.pot.group
      ) %>% factor()
    )  

clean_census_data = clean_census_data %>%
  mutate(
      treat_dist = paste0(
      "treat: ", 
      assigned.treatment,
      ", dist: ", dist.pot.group
      ) %>% factor()
  )  
  



social_perception_baseline = baseline.data %>% 
  select(assigned.treatment, dist.pot.group, cluster.id, county, matches("^(praise|stigma)_[^_]+$")) %>% 
  gather(key = key, value = response, -assigned.treatment, -dist.pot.group, -cluster.id, -county)  %>%
  separate(key, c("praise.stigma", "topic"), "_") %>% 
  separate(topic, c("topic", "question.group"), -2) %>% 
  filter(!is.na(response))  %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  mutate(response_yes = response == "yes") 

baseline.data %>% 
  select(assigned.treatment, dist.pot.group, cluster.id, county, matches("^(praise|stigma)_[^_]+$")) %>%
  select(praise_immunizeA, praise_immunizeB)

praise_df = social_perception_baseline %>%
  filter(topic %in% c("immuniz", "dewor")) %>%
  pivot_wider(
     values_from = response_yes,
     names_from = c(praise.stigma, topic)
  )  %>%
  unnest(c(praise_immuniz, praise_dewor)) 

stigma_df = social_perception_baseline %>%
  filter(topic %in% c("immuniz", "dewor")) %>%
  pivot_wider(
     values_from = response_yes,
     names_from = c(praise.stigma, topic)
  )   %>%
  unnest(c(stigma_immuniz, stigma_dewor))

#### Endline
endline_worm_data = endline.data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  clean_worm_covariates()

endline_implementation_data = endline.data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  clean_implementation_vars()

endline_and_baseline_worm_data = bind_rows(
  endline_worm_data %>%
    select(
      treat_dist, 
      any_of(worm_vars),
      cluster.id, 
      county,
      fully_aware_externalities,
      ) %>%
    mutate(
      type = "endline"
    ),
  baseline_worm_data %>%
    select(
      treat_dist, 
      any_of(worm_vars),
      cluster.id,
      county,
      fully_aware_externalities
      ) %>%
    mutate(
      type = "baseline"
    )
)

endline_vars = c(
  "fully_aware_externalities",
  "all_can_get_worms", 
  "correct_when_treat", 
  "know_worms_infectious"
  )

## If at cluster level, aggregate
if (script_options$community_level) {
  baseline_worm_data = baseline_worm_data %>%
    group_by(cluster.id) %>%
    summarise(
      across(c(worm_vars, cluster.dist.to.pot), mean, na.rm = TRUE),
      treat_dist = unique(treat_dist),
      county = unique(county)
      )


  endline_worm_data = endline_worm_data %>%
    group_by(cluster.id) %>%
    summarise(
      across(endline_vars, mean, na.rm = TRUE),
      treat_dist = unique(treat_dist),
      county = unique(county)
      )

  endline_implementation_data = endline_implementation_data %>%
    group_by(cluster.id) %>%
    summarise(
      across(implementation_vars, mean, na.rm = TRUE),
      treat_dist = unique(treat_dist),
      county = unique(county)
      )

  clean_census_data = clean_census_data %>%
    group_by(cluster.id) %>%
    summarise(
      across(any_of(unique(census_vars)), mean, na.rm = TRUE),
      treat_dist = unique(treat_dist),
      county = unique(county)
    ) 

  analysis_data = analysis_data %>%
    group_by(cluster.id) %>%
    summarise(
      across(any_of(unique(takeup_vars)), mean, na.rm = TRUE),
      treat_dist = unique(treat_dist),
      county = unique(county)
      )
}
## Fits
endline_worm_fit = feols(
    data = endline_worm_data, 
    .[endline_vars] ~ 0 + treat_dist + i(county, ref = "Busia"), 
    cluster = if(script_options$community_level) NULL else ~cluster.id,
    vcov = if(script_options$community_level) "hetero"
    ) 

endline_implementation_fit = feols(
    data = endline_implementation_data, 
    .[implementation_vars] ~ 0 + treat_dist + i(county, ref = "Busia"), 
    cluster = if(script_options$community_level) NULL else ~cluster.id,
    vcov = if(script_options$community_level) "hetero"
    ) 

baseline_worm_fit = feols(
    data = baseline_worm_data, 
    .[worm_vars] ~ 0 + treat_dist + i(county, ref = "Busia"), 
    cluster = if(script_options$community_level) NULL else ~cluster.id,
    vcov = if(script_options$community_level) "hetero"
    ) 

pretreat_fit = feols(
  data = pretreat_data,
  .[pretreat_vars[pretreat_vars != 'religion_christianity']] ~ 0 + treat_dist + i(county, ref = "Busia"),
  cluster = if(script_options$community_level) NULL else ~cluster.id,
  vcov = if(script_options$community_level) "hetero"
)

sens_fit = feols(
  data = clean_sens_summ_imp_df,
  .[sens_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
  vcov =  "hetero"
)




praise_fit = feols(
  praise_df,
  .[praise_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
  cluster = ~cluster.id
)

stigma_fit = feols(
  stigma_df,
  .[stigma_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
  cluster = ~cluster.id
)

# control mean >1 if we have county FE due to LPM, so just don't use county FE.
# This is stupid but the world we live in.
pretreat_christ_fit = feols(
  data = pretreat_data,
  religion_christianity ~ 0 + treat_dist,
  cluster = if(script_options$community_level) NULL else ~cluster.id,
  vcov = if(script_options$community_level) "hetero"
)

census_fit = feols(
    data = clean_census_data, 
    .[census_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
    cluster = if(script_options$community_level) NULL else ~cluster.id,
    vcov = if(script_options$community_level) "hetero"
    ) 

misc_fit = feols(
    data = analysis_data %>%
      select(any_of(takeup_vars), treat_dist, county, cluster.id), 
    .[takeup_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
    cluster = if(script_options$community_level) NULL else ~cluster.id,
    vcov = if(script_options$community_level) "hetero"
  )
# put all the baseline balance fits into a list we can map over
balance_fits = c(
  baseline_worm_fit,
  pretreat_fit,
  list("lhs: religion_christianity" = pretreat_christ_fit),
  census_fit,
  endline_implementation_fit,
  misc_fit,
  sens_fit,
  praise_fit,
  stigma_fit
)


create_balance_comparisons = function(fit) {
  comp_df = avg_comparisons(
    fit,
    variables = list("treat_dist" = "all")
  ) %>%
  as_tibble()


  comp_df = comp_df %>%
    mutate(
      lhs_treatment = str_extract(contrast, "(?<=^treat: )\\w+"), 
      rhs_treatment = str_extract(contrast, "(?<=- treat: )\\w+"), 
      lhs_dist = str_extract(contrast, "(?<=dist: )\\w+"),
      rhs_dist = str_extract(contrast, "(?<=, dist: )\\w+$")
    ) 


  same_dist_subset_comp_df = comp_df  %>%
    filter(
      lhs_dist == rhs_dist,
      rhs_treatment == "control" | lhs_treatment == "control",
      lhs_treatment != rhs_treatment
    ) 

  same_dist_bra_cal_comp_df = comp_df %>%
    filter(lhs_dist == rhs_dist) %>%
    filter(str_detect(contrast, "bracelet") & str_detect(contrast, "calendar"))

    rhs_control_comp_df = same_dist_subset_comp_df %>%
      filter(rhs_treatment == "control") %>%
      bind_rows(
        same_dist_bra_cal_comp_df %>%
          filter(rhs_treatment == "calendar")
      )

    lhs_control_comp_df = same_dist_subset_comp_df %>%
      filter(rhs_treatment != "control") %>%
      bind_rows(
        same_dist_bra_cal_comp_df %>%
          filter(rhs_treatment == "bracelet")
      )

    lhs_control_comp_df = lhs_control_comp_df %>%
      mutate(
        new_estimate = estimate*-1, 
        new_statistic = statistic*-1, 
        new_conf.low = conf.high*-1,
        new_conf.high = conf.low*-1,
        new_lhs_treatment = rhs_treatment,
        new_rhs_treatment = lhs_treatment
      )  %>%
      mutate(
        estimate = new_estimate, 
        statistic = new_statistic,
        conf.low = new_conf.low,
        conf.high = new_conf.high, 
        lhs_treatment = new_lhs_treatment,
        rhs_treatment = new_rhs_treatment
      ) %>%
      select(-contains('new_'))

    rearranged_comp_df = bind_rows(
      lhs_control_comp_df, 
      rhs_control_comp_df
   ) %>%
   select(-contrast)


    control_mean_df = fit %>%
      tidy(conf.int = TRUE) %>%
      filter(str_detect(term, "control")) %>%
      mutate(
        lhs_treatment = "control", rhs_treatment = NA, 
        lhs_dist = if_else(str_detect(term, "close"), "close", "far"), 
        rhs_dist = lhs_dist
      ) %>%
      select(
        -term
      )

  rearranged_comp_df = rearranged_comp_df %>%
    bind_rows(
      control_mean_df
    ) %>%
    mutate(comp_type = "treatment")

    ## Now within treatment across distances

    dist_control_mean_df = fit %>%
      tidy(conf.int = TRUE) %>%
      filter(str_detect(term, "close")) %>%
      mutate(
        lhs_treatment = str_extract(
          term, 
          "(?<=treat: )\\w+"),
        rhs_treatment = NA,
        lhs_dist = "close",
        rhs_dist = NA
      )

    dist_comp_df = comp_df %>%
      filter(
        rhs_dist != lhs_dist,
        lhs_treatment == rhs_treatment
      )  %>%
      select(-contrast) %>%
      bind_rows(
        dist_control_mean_df
      ) %>%
      mutate(comp_type = "distance")

    final_clean_comp_df = bind_rows(
      rearranged_comp_df,
      dist_comp_df
    )
    return(final_clean_comp_df)
}



comp_balance_tidy_df = balance_fits %>%
  map_dfr(
    create_balance_comparisons, 
    .id = "lhs"
  ) 


balance_tidy_df = balance_fits %>%
    map_dfr(tidy, .id = "lhs") %>%
    select(
        lhs, term, estimate, std.error, p.value
    )  

comp_balance_tidy_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "comp_balance_tidy_df.csv"
        )
    )

balance_tidy_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "balance_tidy_df.csv"
        )
    )

construct_joint_test_m = function(object) {
  n_coef = length(coef(object))
  diag_m = diag(n_coef - 1)
  neg_1_m = matrix(-1, nrow = n_coef - 1, ncol = 1)

  hyp_m = cbind(neg_1_m, diag_m)
  return(hyp_m)
}



#### Joint Tests ####
#| joint-tests

## We want to test for balance across all conditions and balance within distance condition
## I don't know how to do such a joint test in R easily so we setup the test matrix 
## manually for the wald test
# Number of dist groups x treatment
n_variables = 8
# matrix R for test 
hyp_matrix = cbind(
  matrix(-1, nrow = n_variables - 1, ncol = 1 ), 
  diag(x = 1, nrow = n_variables - 1, ncol = n_variables)[, 1:(n_variables - 1)]
)

zero_matrix = matrix(0, nrow = 3, ncol = n_variables - 1) 
part_hyp_matrix = zero_matrix
for (i in 1:3) {
  part_hyp_matrix[i, 2*i] = 1
}

hyp_matrix_close = cbind(
  matrix(-1, nrow = 3, ncol = 1), 
  part_hyp_matrix
)

hyp_matrix_far = cbind(
  matrix(0, nrow = 3, ncol = 1),
  matrix(-1, nrow = 3, ncol = 1), 
  part_hyp_matrix[, 1:(ncol(part_hyp_matrix) - 1)]
)

perform_balance_joint_test = function(fit, var, joint_R, close_R, far_R) {
  county_0_mat = matrix(
    0,
    nrow = max(nrow(joint_R), nrow(close_R), nrow(far_R)),
    ncol = coef(fit) %>% length() - 8
    )

  resid_df = fixest::degrees_freedom(fit, type = "resid")
  close_test = car::lht(
    fit,
    cbind(close_R, county_0_mat[1:nrow(close_R), ]),
    error.df = resid_df,
    test = "F"
  )

  far_test = car::lht(
    fit,
    cbind(far_R, county_0_mat[1:nrow(far_R), ]),
    error.df = resid_df,
    test = "F"
  )

  joint_test = car::lht(
    fit,
    cbind(joint_R, county_0_mat[1:nrow(joint_R),]),
    error.df = resid_df,
    test = "F"
  )


  pvals = lst(
    joint_pval = joint_test$`Pr(>F)`[2],
    far_pval = far_test$`Pr(>F)`[2],
    close_pval = close_test$`Pr(>F)`[2]
  ) 

  return(pvals)
}

balance_joint_tests = map(
  balance_fits,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)



balance_joint_tests %>%
    saveRDS(
        file.path(
            script_options$output_path,
            "cluster_balance_fits.rds"
        )
    )

#| baseline-learning

# Here we test if there's a difference in learning between baseline and endline 
# across a range of alternative mechanisms.
comp_endline_vars = endline_vars %>%
  str_remove(., "endline_")
# comp_endline_vars = comp_endline_vars[comp_endline_vars != "know_deworming_stops_worms"]

endline_and_baseline_worm_data %>%
  summarise(across(comp_endline_vars, var))

endline_and_baseline_worm_data %>%
  count(type)

baseline_endline_worm_fit = feols(
      data = endline_and_baseline_worm_data, 
      .[comp_endline_vars] ~ 0 + treat_dist:type + i(county, ref = "Busia"), 
      ~cluster.id
      ) 

perform_worm_change_test = function(data, var, R) {
  subset_data = data %>%
    select(all_of(var), treat_dist, type, cluster.id, county) %>%
      na.omit()

  fml = as.formula(paste0(var, " ~ 0 + treat_dist:type + i(county, ref = 'Busia')"))

  fit = feols(
    data = subset_data,
    fml = fml,
    ~cluster.id
  )
  # n_county = sum(str_detect(names(coef(fit)), "county"))
  nrow_mat = ifelse(is.matrix(R), nrow(R), 1)
  # county_0_mat = matrix(
  #   0,
  #   nrow = nrow_mat,
  #   ncol = n_county
  #   )

  R = matrix(R, nrow = nrow_mat)
  resid_df = fixest::degrees_freedom(fit, type = "resid")
  test = car::lht(
    fit,
    R,
    error.df = resid_df,
    test = "F"
  )
  test_pval = test$`Pr(>F)`[2]
  return(test_pval)
}



indiv_worm_comp_test = baseline_endline_worm_fit %>%
  map(
    ~comparisons(
      .x,
      variable = list(type = "reference"), 
      newdata = datagrid(
        treat_dist = unique(cluster_treat_df$treat_dist)
      )
    )
  )

indiv_worm_comp_df = tibble(
  lhs = names(indiv_worm_comp_test), 
  p.value = map(indiv_worm_comp_test, "p.value"), 
  treat_dist = list(unique(cluster_treat_df$treat_dist)) 
) %>%
  mutate(lhs = str_remove(lhs, "lhs: ")) %>%
    unnest(c(p.value, treat_dist))

#' Another way to generate the hypothesis matrix - slightly more general
generate_joint_worm_hyp_m = function(fit, treat_term, dist_term) {
  hyp_df = fit %>%
    tidy() %>%
    select(term) %>% 
    mutate(
      treat = str_extract(term, "(?<=treat: ).*(?=,)"),
      dist = str_extract(term, "(?<=dist: ).*(?=:)"), 
      type = str_extract(term, "(?<=type).*$")
    ) %>%
    mutate(
      val = 0,
      val = if_else(
        treat == treat_term &
        dist == dist_term &
        type == "baseline", 
        -1, 
        val
        ),
      val = if_else(
        treat == treat_term &
        dist == dist_term &
        type == "endline", 
        1, 
        val
        ),
      val = if_else(str_detect(term, "county::"), 0, val)
    )
  return(hyp_df$val)
}

dist_treat_grid = expand_grid(
  treat = c("bracelet", "calendar", "ink", "control"), 
  dist = c("close", "far")
) %>%
  arrange(dist)

worm_joint_hyp_matrix = map2(
  dist_treat_grid$treat, 
  dist_treat_grid$dist, 
  ~generate_joint_worm_hyp_m(fit = baseline_endline_worm_fit[[1]], .x, .y )
) %>%
  do.call(rbind, .)


joint_worm_p_value =  map_dbl(
  comp_endline_vars,
  ~perform_worm_change_test(
    data = endline_and_baseline_worm_data, 
    var = .x, 
    R = worm_joint_hyp_matrix
  )
)

gen_close_p_val = function(x){
  map(
      c(split(worm_joint_hyp_matrix[1:4, ], 1:4), list(worm_joint_hyp_matrix[1:4, ])),
      ~perform_worm_change_test(endline_and_baseline_worm_data, x, .x)
  )
}


indiv_close_worm_p_value =  map(
  comp_endline_vars, 
  gen_close_p_val
)

gen_far_p_val = function(x) {
  map(
    c(split(worm_joint_hyp_matrix[5:8, ], 1:4), list(worm_joint_hyp_matrix[5:8, ])),
    ~perform_worm_change_test(endline_and_baseline_worm_data,x,  .x)
)
}

indiv_far_worm_p_value = map(comp_endline_vars, gen_far_p_val)


endline_p_val_df = map(
  1:length(comp_endline_vars), 
  ~c(comp_endline_vars[.x], indiv_close_worm_p_value[[.x]] %>% unlist(), indiv_far_worm_p_value[[.x]] %>% unlist(), joint_worm_p_value[[.x]]))  %>%
map_dfr(~data.frame(t(.x)))

treat_levels_c = c("control", "ink", "calendar", "bracelet")
treat_levels = c("ink", "calendar", "bracelet")
col_order = c(
  "lhs", 
  paste0(treat_levels_c, "_close"),
  "close_joint_p",
  paste0(treat_levels_c, "_far"),
  "far_joint_p",
  "joint_p"
)

colnames(endline_p_val_df) = col_order

endline_p_val_df = endline_p_val_df %>%
  mutate(across(c(everything(), -lhs), as.numeric)) %>%
  mutate(fit_type = "pval") 


endline_tidy_df = endline_worm_fit %>%
    map_dfr(tidy, .id = "lhs") %>%
    select(
        lhs, term, estimate, std.error, p.value
    )  


endline_tidy_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "endline_balance_tidy_df.csv"
        )
    )

endline_p_val_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "endline_balance_p_val_df.csv"
        )
    )


#### Social Perception Balanced ####
praise_baseline_fit = social_perception_baseline %>%
  filter(praise.stigma == "praise") %>%
  feols(
    response_yes ~ 0 + treat_dist + i(county, ref = "Busia"),
    ~cluster.id,
    split = ~topic
  )

stigma_baseline_fit = social_perception_baseline %>%
  filter(praise.stigma == "stigma") %>%
  feols(
    response_yes ~ 0 + treat_dist + i(county, ref = "Busia"),
    ~cluster.id,
    split = ~topic
  )

praise_stigma_fits = c(
  praise_baseline_fit,
  stigma_baseline_fit
)

praise_stigma_joint_tests  = map(
  praise_stigma_fits,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

saveRDS(
  list(
    praise_stigma_fits,
    praise_stigma_joint_tests
  ),
  file.path(
      script_options$output_path,
      "praise_stigma_baseline.rds"
  )
)




#### Continuous Distance Tests ####
baseline_worm_dist_fit = feols(
  data = baseline_worm_data,
  .[worm_vars] ~ 0 + cluster.dist.to.pot + i(county, ref = "Busia"), 
  ~cluster.id
)

pretreat_dist_fit = feols(
  data = pretreat_data,
  .[pretreat_vars] ~ 0 + cluster.dist.to.pot + i(county, ref = "Busia"), 
  ~cluster.id 
)

census_dist_fit = feols(
    data = clean_census_data, 
    .[census_vars] ~ 0 + cluster.dist.to.pot + i(county, ref = "Busia"),
    cluster = ~cluster.id
    ) 

ri_fun = function(draw) {
  set.seed(draw)
  perm_baseline_worm_data = baseline_worm_data %>%
    mutate(
      perm_dist = sample(cluster.dist.to.pot, size = n())
    )
  perm_census_data = clean_census_data %>%
    mutate(perm_dist = sample(cluster.dist.to.pot, size = n()))
  perm_pretreat_data = pretreat_data %>%
    mutate(perm_dist = sample(cluster.dist.to.pot, size = n())) 

  perm_baseline_worm_dist_fit = feols(
    data = perm_baseline_worm_data,
    .[worm_vars] ~ 0 + perm_dist + i(county, ref = "Busia"), 
    cluster = ~cluster.id
  )
  perm_census_dist_fit = feols(
      data = perm_census_data, 
      .[census_vars] ~ 0 + perm_dist + i(county, ref = "Busia"),
      cluster = ~cluster.id
      ) 
  perm_pretreat_dist_fit = feols(
    data = perm_pretreat_data,
    .[pretreat_vars] ~ 0 + perm_dist + i(county, ref = "Busia"), 
    ~cluster.id
  )

  
  perm_coef_df =  
    bind_rows(
      map_dfr(perm_baseline_worm_dist_fit, tidy, .id = "lhs") %>%
        filter(term == "perm_dist"),
      map_dfr(perm_census_dist_fit, tidy, .id = "lhs") %>%
        filter(term == "perm_dist"),
      map_dfr(perm_pretreat_dist_fit, tidy, .id = "lhs") %>%
        filter(term == "perm_dist")
    ) %>% 
      mutate(draw = draw)

  return(
    perm_coef_df
  )
}

if (script_options$fit_ri) {
  plan(multisession, workers = 12)
  perm_fit_df = future_map_dfr(
    1:100, 
    ri_fun, 
    .progress = TRUE, 
    .options = furrr_options(
      seed = TRUE,
      packages = c("broom", "fixest")
      )
    )
  saveRDS(perm_fit_df, "temp-data/balance-cts-dist-ri.rds")
} else {
  perm_fit_df = read_rds("temp-data/balance-cts-dist-ri.rds")
}

lhs_translation_df = tribble(
  ~lhs, ~clean_name,
  "lhs: age", "Age",
  "lhs: all_can_get_worms", "Know everyone can be infected",
  "lhs: treated_lgl", "Dewormed in the past",
  "lhs: floor_tile_cement", "Floor made of tile/cement",
  "lhs: completed_primary", "Completed primary schooling",
  "lhs: correct_when_treat", "Know bi-yearly treatment recommended",
  "lhs: fully_aware_externalities", "Understands externalities",
  "lhs: female", "Female",
  "lhs: phone_owner", "Phone owner",
  "lhs: know_deworming_stops_worms", "Knows deworming stops worms"
)
realised_fit_df = bind_rows(
  map_dfr(baseline_worm_dist_fit, tidy, .id = "lhs"),
  map_dfr(census_dist_fit, tidy, .id = "lhs"),
  map_dfr(pretreat_dist_fit, tidy, .id = "lhs")
) %>%
  filter(term == "cluster.dist.to.pot")

realised_fit_df = realised_fit_df %>%
  left_join(
    lhs_translation_df,
    by = "lhs"
  )

realised_fit_df %>%
  select(lhs, statistic)

plot_perm_fit_df = perm_fit_df %>%
  left_join(
    realised_fit_df %>%
      select(lhs, realised_statistic = statistic), 
    by = "lhs" )  %>%
  left_join(
    lhs_translation_df,
    by = "lhs"
  )


ri_p_val_df = plot_perm_fit_df %>%
  group_by(clean_name) %>%
  summarise(
    p_val = paste0("p = ", round(mean(statistic > realised_statistic), 3)),
    realised_statistic = unique(realised_statistic),
    x = quantile(statistic, 0.95, na.rm = TRUE)
  ) 


balance_data = lst(
  analysis_data, 
  pretreat_data, 
  baseline_worm_data,
  endline_worm_data,
  endline_and_baseline_worm_data,
  clean_census_data,
  pretreat_data,
  endline_vars, 
  worm_vars,
  pretreat_vars,
  census_vars,
  plot_perm_fit_df,
  ri_p_val_df
)


comp_balance_tidy_df %>%
  filter(term == "treat_dist") %>%
  select(
    lhs, 
    lhs_treatment, 
    rhs_treatment, 
    lhs_dist, 
    rhs_dist, 
    comp_type,
    estimate,
    p.value
    )  %>%
    write_csv(
      file.path(
        script_options$output_path,
        "comm-pvals.csv"
      )
    )


tibble(
  lhs = names(balance_joint_tests),
  joint = map_dbl(balance_joint_tests, "joint_pval"),
  close = map_dbl(balance_joint_tests, "close_pval"),
  far = map_dbl(balance_joint_tests, "far_pval")
) %>%
  mutate(across(where(is.numeric), round, 5))


saveRDS(
  balance_data, 
  file.path(
    script_options$output_path,
    "saved_balance_data.rds"
  )
)




# #### Continuous Balance - regressions
# #| cts-dist-balance
# dist_balance_cts_fun = function(rhs_var) {
#   ## Indiv
#   indiv_balance_cts_fit = feols(
#       data = clean_census_data %>%
#         mutate(dist_measure = {{ rhs_var }}/1000), 
#       .[indiv_balance_vars] ~  dist_measure + dist_measure^2,
#       cluster = ~cluster.id
#       ) 

#   ## Know

#   ## baseline
#   baseline_balance_cts_fit = feols(
#       data = baseline_balance_data %>%
#         mutate(dist_measure = {{ rhs_var }}/1000), 
#       .[baseline_vars] ~  dist_measure + dist_measure^2, 
#       ~cluster.id
#       ) 

#     return(
#       c(
#         indiv_balance_cts_fit,
#         baseline_balance_cts_fit
#       )
#     )
# }

# dist_balance_disc_fun = function(rhs_var, interval_length) {
#   clean_census_data = clean_census_data %>%
#     mutate(
#       dist_measure = cut_interval({{ rhs_var }}, length = interval_length )
#     )

#   baseline_balance_data = baseline_balance_data %>%
#     mutate(
#       dist_measure = cut_interval({{ rhs_var }}, length = interval_length )
#     )

#   ## Indiv
#   indiv_balance_disc_fit = feols(
#       data = clean_census_data,
#       .[indiv_balance_vars] ~  0 + dist_measure + i(county, ref = "Busia"),
#       cluster = ~cluster.id
#       ) 

#   ## baseline
#   baseline_balance_disc_fit = feols(
#       data = baseline_balance_data,
#       .[baseline_vars] ~ 0 + dist_measure + i(county, ref = "Busia"), 
#       ~cluster.id
#       ) 

#     return(
#       c(
#         indiv_balance_disc_fit,
#         baseline_balance_disc_fit
#       )
#     )
# }


# disc_dist_balance = dist_balance_disc_fun(cluster.dist.to.pot, interval_length = 750)
# fully_cts_dist_balance = dist_balance_cts_fun(cluster.dist.to.pot)


# fully_cts_dist_balance %>%
#   saveRDS(
#     file.path(
#       script_options$output_path,
#       "fully_cts_dist_balance.rds"
#     )
#   )


# construct_joint_test_m = function(object) {
#   n_coef = length(coef(object)) - 2
#   diag_m = diag(n_coef - 1)
#   neg_1_m = matrix(-1, nrow = n_coef - 1, ncol = 1)
#   mat_0 = matrix(0, nrow = n_coef -1, ncol = 2)
#   hyp_m = cbind(neg_1_m, diag_m)
#   hyp_m = cbind(hyp_m, mat_0)
#   return(hyp_m)
# }


# R_dist_cts_discrete = construct_joint_test_m(disc_dist_balance[[1]])
# pval_cts_binned_dist = map_dbl(
#   disc_dist_balance,
#   ~car::lht(
#     .x,
#     R_dist_cts_discrete,
#     test = "F"
#   )$`Pr(>Chisq)`[2]
# ) 

# pval_cts_binned_dist %>%
#   saveRDS(
#     file.path(
#       script_options$output_path,
#       "discrete_cts_dist_balance.rds"
#     )
#   )

