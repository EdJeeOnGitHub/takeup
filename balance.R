#!/usr/bin/Rscript
script_options = docopt::docopt(
  stringr::str_glue("Usage:
  balance.R [options]
  
Options:
  --output-path=<path>  Where to save output files [default: {file.path('temp-data')}]
  --community-level
  --fit-ri
  --orig
  --attrition
  --monitored-attrition
  --sms
  --main

"),
  args = if (interactive()) "
    --output-path=temp-data \
    --sms
    " else commandArgs(trailingOnly = TRUE)
)


run_all <- !any(
  script_options$main, 
  script_options$attrition, 
  script_options$monitored_attrition, 
  script_options$continuous_distance, 
  script_options$orig, 
  script_options$fit_ri,
  script_options$sms
)

set.seed(12932)

construct_joint_test_m = function(object) {
  n_coef = length(coef(object))
  diag_m = diag(n_coef - 1)
  neg_1_m = matrix(-1, nrow = n_coef - 1, ncol = 1)

  hyp_m = cbind(neg_1_m, diag_m)
  return(hyp_m)
}




library(tidyverse)
library(marginaleffects)
library(broom)
library(knitr)
library(kableExtra)
library(ggthemes)
library(fixest)
library(magrittr)
library(furrr)

source(file.path("scratch", "reduced-form-setup.R"))
source("balance-functions.R")
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


lhs_translation_df = tribble(
  ~lhs, ~clean_name, ~sample,
  "lhs: age.census", "Age", "takeup_sample",
  # "lhs: all_can_get_worms", "Know everyone can be infected", "baseline_worm",

  "lhs: adults_can_get_worms", "Know adults get worms", "baseline_worm",
  "lhs: know_children_get_worms", "Know children get worms", "baseline_worm",
  "lhs: sick_worms_only", "Believe deworming is for the sick", "baseline_worm",
  "lhs: know_medicine_stops_worms", "Know medication treats worms", "baseline_worm",
  "lhs: treated_lgl", "Dewormed in the past", "baseline_worm",
  "lhs: treated_past_year", "Dewormed in the past year", "baseline_worm",
  "lhs: correct_when_treat", "Know bi-yearly treatment recommended", "baseline_worm",
  "lhs: externality_omnibus", "Know worms impose externality", "baseline_worm",


  "lhs: floor_tile_cement", "Floor made of tile/cement", "pretreat",
  "lhs: completed_primary", "Completed primary schooling", "pretreat",
  "lhs: female", "Female", "takeup_sample",
  "lhs: have_phone_lgl", "Phone owner", "takeup_sample",
  "lhs: n_per_cluster", "Number of individuals per community", "takeup_sample",
  "lhs: cluster.dist.to.pot", "Distance to PoT", "takeup_sample",
  "lhs: years_schooling", "Years schooling", "pretreat",
  "lhs: ethnicity_luhya", "Main ethnicity/Luhya", "pretreat",
  "lhs: religion_christianity", "Christian", "pretreat",
  "lhs: know_deworm", "Know about community-based MDA?", "implementation",
  "lhs: treat_begin", "Know when MDA starts?", "implementation",
  "lhs: treat_end", "Know about MDA ends?", "implementation",
  "lhs: days_available", "Know length of MDA?", "implementation",
  "lhs: chv_visit", "Did a CHV visit you?", "implementation",

  "lhs: stigma_dewor", "Would you judge: Not deworming?", "social_image_concerns",
  "lhs: stigma_immuniz", "Would you judge: Not immunizing a child?", "social_image_concerns",
  "lhs: praise_dewor", "Would you praise: Deworming?", "social_image_concerns",
  "lhs: praise_immuniz", "Would you praise: Immunizing a child?", "social_image_concerns",

  "lhs: pct_announce", "Announcement about MDA in your community?", "implementation",
  "lhs: pct_knowledge_message", "CHV share deworming practices?", "implementation",
  "lhs: pct_avail_message", "CHV share where to get dewormed?", "implementation",


  "lhs: belief_sample_n_per_cluster", "Number of individuals per community", "endline_belief_sample",
  "lhs: belief_sample_female", "Female", "endline_belief_sample",
  "lhs: belief_sample_have_phone_lgl", "Phone owner", "endline_belief_sample",
  "lhs: belief_sample_age.census", "Age", "endline_belief_sample",
  "lhs: belief_sample_cluster.dist.to.pot", "Distance to PoT", "endline_belief_sample"
) 


worm_vars = c(
  "adults_can_get_worms", 
  "know_children_get_worms",
  "sick_worms_only",
  "know_medicine_stops_worms",
  "treated_lgl",
  "treated_past_year",
  "correct_when_treat",
  "externality_omnibus"
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
  mutate(cluster.id = factor(cluster.id)) %>%
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
  
social_perception_baseline = baseline_data %>% 
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
baseline_worm_data = baseline_data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  clean_worm_covariates()

endline_worm_data = endline_data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  clean_worm_covariates()

endline_implementation_data = endline_data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  clean_implementation_vars()

endline_vars = c(
  "fully_aware_externalities",
  "all_can_get_worms", 
  "correct_when_treat", 
  "know_worms_infectious"
  )

# Removing list columns as they mess up marginaleffects
for (nm in c("endline_worm_data", "endline_implementation_data", "baseline_worm_data",
             "pretreat_data", "clean_sens_imp_df", "praise_df", "stigma_df",
             "clean_census_data", "analysis_data")) {
  assign(nm, select(get(nm), where(~!is.list(.))))
}

summ_know_A_df = summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A"))

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
# Get the keys for the belief sample so we can do balance tests on this sample
endline_belief_keys = endline_know_A_df$KEY.individ

clust_n_df = n_indiv_df %>%
  select(cluster.id, n_per_cluster)

endline_know_balance_df = census_data %>%
  filter(KEY.individ %in% endline_belief_keys) %>%
  left_join(
    clust_n_df,
    by = "cluster.id"
  ) %>%
  left_join(
    cluster_treat_df,
    by = "cluster.id"
  ) %>%
  mutate(
    female = gender == 2,
    age = age.census,
    have_phone_lgl = have_phone == "Yes" 
  ) 



if (run_all || script_options$main) {

takeup_balance_input = analysis_data %>%
  select(any_of(takeup_vars), treat_dist, county, cluster.id) %>%
  # convert to Km for this table
  mutate(cluster.dist.to.pot = cluster.dist.to.pot / 1000 )

endline_know_balance_input = endline_know_balance_df %>%
  select(any_of(takeup_vars), treat_dist, county, cluster.id) %>%
  # convert to Km for this table
  mutate(cluster.dist.to.pot = cluster.dist.to.pot / 1000 )

takeup_bal_fit = feols(
    data = takeup_balance_input,
    .[takeup_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
    cluster = ~cluster.id
  )

endline_know_bal_fit = feols(
  data = endline_know_balance_input,
    .[takeup_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
    cluster = ~cluster.id
)


takeup_bal_joint_test = map(
  takeup_bal_fit,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

endline_know_bal_joint_test = map(
  endline_know_bal_fit,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

simple_balance = list(
  takeup_bal_fit = takeup_bal_fit,
  endline_know_bal_fit = endline_know_bal_fit,
  takeup_bal_joint_test = takeup_bal_joint_test,
  endline_know_bal_joint_test = endline_know_bal_joint_test
)

saveRDS(
  simple_balance,
  file.path(
    script_options$output_path,
    "main-balance-samples.rds"
  )
)


comp_takeup_balance_tidy_df = takeup_bal_fit %>%
  map_dfr(
    ~create_balance_comparisons(.x, data = takeup_balance_input),
    .id = "lhs",
    .progress = TRUE
  )
comp_endline_know_balance_tidy_df = endline_know_bal_fit %>%
  map_dfr(
    ~create_balance_comparisons(.x, data = endline_know_balance_input),
    .id = "lhs",
    .progress = TRUE
  )

write_csv(
  comp_takeup_balance_tidy_df,
  file.path(
    script_options$output_path,
    "comp_takeup_balance_tidy_df.csv"
  )
)

write_csv(
  comp_endline_know_balance_tidy_df,
  file.path(
    script_options$output_path,
    "comp_endline_know_balance_tidy_df.csv"
  )
)

takeup_bal_fit %>%
  map_dfr(tidy, .id = "lhs") %>%
  select(
    lhs, term, estimate, std.error, p.value
  ) %>%
  write_csv(
    file.path(
      script_options$output_path,
      "takeup_balance_tidy_df.csv"
    )
  )

endline_know_bal_fit %>%
  map_dfr(tidy, .id = "lhs") %>%
  select(
    lhs, term, estimate, std.error, p.value
  ) %>%
  write_csv(
    file.path(
      script_options$output_path,
      "endline_know_balance_tidy_df.csv"
    )
  )
}

if (run_all || script_options$orig)  {
## If at cluster level, aggregate
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
      select(any_of(takeup_vars), treat_dist, county, cluster.id) %>%
      # convert to Km for this table
      mutate(cluster.dist.to.pot = cluster.dist.to.pot / 1000 ), 
    .[takeup_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
    cluster = if(script_options$community_level) NULL else ~cluster.id,
    vcov = if(script_options$community_level) "hetero"
  )



endline_misc_vars = paste0("belief_sample_", takeup_vars)

misc_fit_endline_belief_sample = feols(
    data = analysis_data %>%
      filter(KEY.individ %in% endline_belief_keys) %>%
      select(any_of(takeup_vars), treat_dist, county, cluster.id) %>%
      # convert to Km for this table
      mutate(cluster.dist.to.pot = cluster.dist.to.pot / 1000 ) %>%
      rename_with(
        ~paste0("belief_sample_", .), any_of(takeup_vars)
      ),
    .[endline_misc_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
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
  stigma_fit,
  misc_fit_endline_belief_sample
)

balance_fit_data = c(
  rep(list(baseline_worm_data), length(baseline_worm_fit)),
  rep(list(pretreat_data), length(pretreat_fit)),
  list(pretreat_data),
  rep(list(clean_census_data), length(census_fit)),
  rep(list(endline_implementation_data), length(endline_implementation_fit)),
  rep(list(analysis_data %>%
      select(any_of(takeup_vars), treat_dist, county, cluster.id) %>%
      mutate(cluster.dist.to.pot = cluster.dist.to.pot / 1000)), length(misc_fit)),
  rep(list(clean_sens_summ_imp_df), length(sens_fit)),
  rep(list(praise_df), length(praise_fit)),
  rep(list(stigma_df), length(stigma_fit)),
  rep(list(analysis_data %>%
      filter(KEY.individ %in% endline_belief_keys) %>%
      select(any_of(takeup_vars), treat_dist, county, cluster.id) %>%
      mutate(cluster.dist.to.pot = cluster.dist.to.pot / 1000) %>%
      rename_with(
        ~paste0("belief_sample_", .), any_of(takeup_vars)
      )), length(misc_fit_endline_belief_sample))
)


comp_balance_tidy_df = balance_fits %>%
  map2_dfr(
    balance_fit_data,
    ~create_balance_comparisons(.x, data = .y),
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



#### Overall Means for Intro ####
baseline_worm %>%
  summarise(
    across(
      c(
        "know_medicine_stops_worms",
        "correct_when_treat",
        "know_children_get_worms",
        "all_can_get_worms",
        "sick_worms_only"
      ),
      mean, na.rm = TRUE
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "mean"
  ) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  write_csv(
    file.path(
      script_options$output_path,
      "baseline_worm_means.csv"
    )
  )

#### Joint Tests ####
#| joint-tests

## We want to test for balance across all conditions and balance within distance condition
## I don't know how to do such a joint test in R easily so we setup the test matrix 
## manually for the wald test
# Number of dist groups x treatment

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


dist_treat_grid = expand_grid(
  treat = c("bracelet", "calendar", "ink", "control"), 
  dist = c("close", "far")
) %>%
  arrange(dist)


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
  .[worm_vars] ~  cluster.dist.to.pot | county, 
  ~cluster.id
)

pretreat_dist_fit = feols(
  data = pretreat_data,
  .[pretreat_vars] ~  cluster.dist.to.pot | county, 
  ~cluster.id 
)

census_dist_fit = feols(
    data = clean_census_data, 
    .[census_vars] ~  cluster.dist.to.pot | county,
    cluster = ~cluster.id
    ) 

takeup_dist_fit = feols(
  data = analysis_data %>%
    select(any_of(takeup_vars), treat_dist, county, cluster.id), 
  .[takeup_vars] ~  cluster.dist.to.pot | county,
  cluster = ~cluster.id
)





if (script_options$fit_ri) {

  load_counterfactual_distance_pools = function() {
    counterfactual_distance_df = read_csv(
      file.path("data", "simulated-counterfactual-treatment-assignment.csv"),
      col_types = cols(
        .default = col_skip(),
        cluster.id = col_character(),
        dist = col_double()
      )
    ) %>%
      filter(!is.na(cluster.id)) %>%
      mutate(
        cluster.id = as.character(as.integer(cluster.id)),
        perm_dist = dist / 1000
      ) %>%
      select(cluster.id, perm_dist)

    split(counterfactual_distance_df$perm_dist, counterfactual_distance_df$cluster.id)
  }

  assign_counterfactual_cluster_dists = function(data, seed, distance_pools) {
    set.seed(seed)
    cluster_ids = data %>%
      mutate(cluster.id = as.character(cluster.id)) %>%
      pull(cluster.id) %>%
      unique()

    missing_clusters = setdiff(cluster_ids, names(distance_pools))
    if (length(missing_clusters) > 0) {
      stop(
        "Missing counterfactual distance simulations for clusters: ",
        paste(missing_clusters, collapse = ", ")
      )
    }

    counterfactual_cluster_dist_df = tibble(
      cluster.id = cluster_ids,
      perm_dist = vapply(
        cluster_ids,
        function(cluster_id) sample(distance_pools[[cluster_id]], 1),
        numeric(1)
      )
    )

      data = data %>%
        mutate(cluster.id = as.character(cluster.id)) %>%
        left_join(
          counterfactual_cluster_dist_df,
          by = "cluster.id"
        )
      return(data)
  }

  counterfactual_distance_pools = load_counterfactual_distance_pools()

  ri_fun = function(draw) {
    set.seed(draw)
    perm_baseline_worm_data = baseline_worm_data %>%
      assign_counterfactual_cluster_dists(draw, counterfactual_distance_pools)
    perm_census_data = clean_census_data %>%
      assign_counterfactual_cluster_dists(draw, counterfactual_distance_pools)

    perm_pretreat_data = pretreat_data %>%
      assign_counterfactual_cluster_dists(draw, counterfactual_distance_pools)

    perm_takeup_data = analysis_data %>%
      select(any_of(takeup_vars), cluster.dist.to.pot, county, cluster.id) %>%
      assign_counterfactual_cluster_dists(draw, counterfactual_distance_pools)

    perm_baseline_worm_dist_fit = feols(
      data = perm_baseline_worm_data,
      .[worm_vars] ~  perm_dist | county, 
      cluster = ~cluster.id
    )
    perm_census_dist_fit = feols(
        data = perm_census_data, 
        .[census_vars] ~  perm_dist | county,
        cluster = ~cluster.id
        ) 
    perm_pretreat_dist_fit = feols(
      data = perm_pretreat_data,
      .[pretreat_vars] ~  perm_dist | county, 
      cluster = ~cluster.id
    )

    perm_takeup_dist_fit = feols(
      data = perm_takeup_data, 
      .[takeup_vars] ~  perm_dist | county,
      cluster = ~cluster.id
    )

    
    perm_coef_df =  
      bind_rows(
        map_dfr(perm_baseline_worm_dist_fit, tidy, .id = "lhs") %>%
          filter(term == "perm_dist"),
        map_dfr(perm_census_dist_fit, tidy, .id = "lhs") %>%
          filter(term == "perm_dist"),
        map_dfr(perm_pretreat_dist_fit, tidy, .id = "lhs") %>%
          filter(term == "perm_dist"),
        map_dfr(
          perm_takeup_dist_fit, 
          tidy, 
          .id = "lhs"
        ) %>%
          filter(term == "perm_dist"
        )
      ) %>% 
        mutate(draw = draw)

    return(
      perm_coef_df
    )
  }

    
  plan(multisession, workers = 12)
  perm_fit_df = future_map_dfr(
    1:500, 
    ri_fun, 
    .progress = TRUE,
    .options = furrr_options(
      seed = TRUE,
      packages = c("broom", "dplyr", "fixest", "tibble")
      )
    )
  saveRDS(perm_fit_df, "temp-data/balance-cts-dist-ri-fe.rds")
} else {
  perm_fit_df = read_rds("temp-data/balance-cts-dist-ri-fe.rds")
}




realised_fit_df = bind_rows(
  map_dfr(baseline_worm_dist_fit, tidy, .id = "lhs"),
  map_dfr(census_dist_fit, tidy, .id = "lhs"),
  map_dfr(pretreat_dist_fit, tidy, .id = "lhs"),
  map_dfr(takeup_dist_fit, tidy, .id = "lhs")
) %>%
  filter(term == "cluster.dist.to.pot")







realised_fit_df = realised_fit_df %>%
  inner_join(
    lhs_translation_df,
    by = "lhs"
  ) 



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
  group_by(lhs) %>%
  summarise(
    p_val = paste0("p = ", round(mean(abs(statistic) > abs(realised_statistic)), 3)),
    realised_statistic = unique(realised_statistic),
    x = quantile(statistic, 0.95, na.rm = TRUE),
    clean_name = unique(clean_name)
  ) 

saveRDS(
  list(
    plot_perm_fit_df = plot_perm_fit_df,
    ri_p_val_df = ri_p_val_df,
    realised_fit_df = realised_fit_df
  ),
  file.path(
    script_options$output_path,
    "balance-ri-fe.rds"
  )
)

N_endline_belief = analysis_data %>%
  filter(KEY.individ %in% endline_belief_keys) %>%
  nrow()

balance_data = lst(
  analysis_data,
  pretreat_data,
  baseline_worm_data,
  endline_worm_data,
  clean_census_data,
  pretreat_data,
  endline_vars,
  worm_vars,
  pretreat_vars,
  census_vars,
  plot_perm_fit_df,
  ri_p_val_df,
  realised_fit_df,
  N_endline_belief
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
}

#### Knowledge Table Attrition Analysis ####
if (run_all || script_options$attrition) {
library(fixest)

endline_data = endline_data %>%
  mutate(
    not_in_know_table = KEY.individ %in% in_endline_not_know_table
  )

endline_data %>%
  count(not_in_know_table, in_know_table)


all_endline_data = all_endline_data %>%
  mutate(
    not_in_know_table = KEY.individ %in% in_endline_not_know_table
  )


nrow(endline_data)
nrow(all_endline_data)
nrow(all_endline_data_frame)

# First attrition table: does treatment predict attrition from the know table?
# (1) Pooled treatment effect on attrition
attrition_pooled = endline_data %>%
  feols(
    not_in_know_table ~ i(assigned.treatment, ref = "control") | county,
    cluster = ~cluster.id
  )
# (2) Treatment x distance interaction
attrition_treat_dist = endline_data %>%
  feols(
    not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") | county,
    cluster = ~cluster.id
  )

# (3) Full sample (include SMS treat/control)
attrition_sms_full = all_endline_data %>%
  feols(
    not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") | county,
    cluster = ~cluster.id
  )

# (3) F-test: joint significance of treatment indicators
attrition_ftest = endline_data %>%
  feols(
    not_in_know_table ~ i(assigned.treatment, ref = "control") | county,
    cluster = ~cluster.id
  ) %>%
  wald("assigned.treatment")

setFixest_dict(c(
  not_in_know_table = "Missing from Know Table",
  "assigned.treatment::ink" = "Ink",
  "assigned.treatment::calendar" = "Calendar",
  "assigned.treatment::bracelet" = "Bracelet",
  gender_num = "Gender",
  have_phone_num = "Has Phone",
  dist.to.pot = "Distance to POT (m)",
  "dist.pot.group" = "Distance Group",
  "close" = "Close",
  "far" = "Far",
  age = "Age"
))

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
# ---- Attrition predicted by treatment status ----
etable(attrition_pooled, attrition_treat_dist, attrition_sms_full,
       headers = c("Pooled", "Treat $\\times$ Distance", "SMS Treatment Included"),
       depvar = FALSE,
       fitstat = c("n", "r2"),
       se.below = TRUE,
       tex = TRUE,
       title = "Differential Attrition from Knowledge Table by Treatment Assignment",
       label = "tab:attrition-treatment",
       notes = "",
       file = "presentations/tables/attrition-by-treatment.tex",
       replace = TRUE,
       postprocess.tex = tex_postprocessing,
       digits = 3,
       digits.stats = 3,
       drop.section = "fixef",
       style.df = style.df(depvar.title = "", fixef.title = "", var.title = "", stats.title = "")
)

# Second attrition table: are attrition housholds different on observables
attrition_endline_data = endline_data %>%
  select(-female, -age.census, -all_of(implementation_vars)) %>%
  left_join(
    select(clean_census_data, KEY.individ, all_of(census_vars)) %>% distinct(),
    by = "KEY.individ"
  ) %>%
  left_join(
    select(analysis_data, KEY.individ, any_of(takeup_vars)) %>% distinct(),
    by = "KEY.individ"
  ) %>%
  left_join(
    select(pretreat_data, KEY.individ, all_of(pretreat_vars)) %>% distinct(),
    by = "KEY.individ"
  ) %>%
  left_join(
    select(endline_worm_data, KEY.individ, any_of(worm_vars)) %>% distinct(),
    by = "KEY.individ"
  ) %>%
  left_join(
    select(endline_implementation_data, KEY.individ, all_of(implementation_vars)) %>% distinct(),
    by = "KEY.individ"
  ) %>%
  left_join(
    select(clean_sens_summ_imp_df, cluster.id, any_of(sens_vars)) %>% distinct()
  )

attrition_balance_vars = c(
  age.census = "Age",
  female = "Female",
  have_phone_lgl = "Phone owner",
  n_per_cluster = "Number of individuals per community",
  cluster.dist.to.pot =  "Distance to PoT", 

  floor_tile_cement = "Floor made of tile/cement",
  completed_primary = "Completed primary schooling",
  ethnicity_luhya = "Main ethnicity/Luhya",
  religion_christianity = "Christian",

  adults_can_get_worms = "Know adults can get worms",
  know_children_get_worms = "Know children get worms",
  sick_worms_only = "Believe deworming is for the sick",
  know_medicine_stops_worms = "Know medication treats worms",
  correct_when_treat = "Know bi-yearly treatment recommended",
  externality_omnibus = "Know worms impose externality",

  chv_visit = "Did a CHV visit you?",
  pct_announce = "Announcement about MDA in your community?",
  pct_knowledge_message = "CHV share deworming practices?",
  pct_avail_message = "CHV share where to get dewormed?"
)


attrition_balance_df = map_dfr(names(attrition_balance_vars), function(v) {
  fit = feols(as.formula(paste0(v, " ~ not_in_know_table | county")),
              data = attrition_endline_data, cluster = ~cluster.id)
  ct = coeftable(fit)

  group_means = attrition_endline_data %>%
    summarise(
      mean_present = mean(.data[[v]][!not_in_know_table], na.rm = TRUE),
      mean_missing = mean(.data[[v]][not_in_know_table], na.rm = TRUE),
      n_present = sum(!is.na(.data[[v]]) & !not_in_know_table),
      n_missing = sum(!is.na(.data[[v]]) & not_in_know_table)
    )

  tibble(
    variable = v,
    label = attrition_balance_vars[v],
    mean_present = group_means$mean_present,
    mean_missing = group_means$mean_missing,
    diff = ct["not_in_know_tableTRUE", "Estimate"],
    se = ct["not_in_know_tableTRUE", "Std. Error"],
    pval = ct["not_in_know_tableTRUE", "Pr(>|t|)"],
    n_present = group_means$n_present,
    n_missing = group_means$n_missing
  )
})

# Console output
attrition_balance_df %>%
  mutate(stars = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**", pval < 0.1 ~ "*", TRUE ~ "")) %>%
  select(label, mean_present, mean_missing, diff, se, pval, stars) %>%
  print(n = Inf)

# ---- .tex output ----
fmt = function(x, d = 3) formatC(x, format = "f", digits = d)
stars_tex = function(p) case_when(p < 0.01 ~ "^{***}", p < 0.05 ~ "^{**}", p < 0.1 ~ "^{*}", TRUE ~ "")

tex_rows = attrition_balance_df %>%
  mutate(row = paste0(
    label, " & ",
    fmt(mean_present), " & ",
    fmt(mean_missing), " & ",
    "$", fmt(diff), stars_tex(pval), "$ & ",
    "$(", fmt(se), ")$ & ",
    format(n_present + n_missing, big.mark = ","),
    " \\\\"
  )) %>%
  pull(row)

n_present_total = max(attrition_balance_df$n_present)
n_missing_total = max(attrition_balance_df$n_missing)

att_demo_vars = c("age.census", "female", "phone_owner", "n_per_cluster", "cluster.dist.to.pot", "floor_tile_cement", "completed_primary", "ethnicity_luhya", "religion_christianity")

att_know_vars = c("adults_can_get_worms", "know_children_get_worms", "sick_worms_only", "know_medicine_stops_worms", "correct_when_treat", "externality_omnibus")

att_impl_vars = c("chv_visit", "pct_announce", "pct_knowledge_message", "pct_avail_message")

tex_table = c(
  "\\centering",
  "\\begin{tabular}{l ccccc}",
  "\\hline\\hline",
  paste0(" & \\multicolumn{1}{c}{Present} & \\multicolumn{1}{c}{Missing} & ",
         "\\multicolumn{1}{c}{Difference} & \\multicolumn{1}{c}{SE} & \\multicolumn{1}{c}{N} \\\\"),
  "\\hline",
  "\\addlinespace",
  "\\textit{Demographics} \\\\",
  "\\addlinespace",
  tex_rows[attrition_balance_df$variable %in% att_demo_vars],
  "\\addlinespace",
  "\\textit{Deworming Knowledge} \\\\",
  "\\addlinespace",
  tex_rows[attrition_balance_df$variable %in% att_know_vars],
  "\\addlinespace",
  "\\textit{Implementation} \\\\",
  "\\addlinespace",
  tex_rows[attrition_balance_df$variable %in% att_impl_vars],
  "\\hline\\hline",
  "\\end{tabular}"
)

writeLines(tex_table, "presentations/tables/attrition-covariate-balance.tex")


# Finally, are baseline characteristics of attrition HHs in treat different to control
# Among attrition HHs: y_vhB = β₀ + β₁ T_vh + ε_vhB, clustered at cluster level

attrition_treat_data = attrition_endline_data %>%
  select(-where(is.list)) %>%
  filter(not_in_know_table) %>%
  mutate(
    treat_dist = paste0("treat: ", assigned.treatment, ", dist: ", dist.pot.group) %>% factor()
  )

n_attrition = nrow(attrition_treat_data)

attrition_treat_fits = feols(
  data = attrition_treat_data,
  .[names(attrition_balance_vars)] ~ 0 + treat_dist + i(county, ref = "Busia"),
  cluster = ~cluster.id
)

# Pairwise comparisons (control means, ink-con, cal-con, bra-con, bra-cal)
attrition_treat_comp_df = attrition_treat_fits %>%
  map_dfr(~create_balance_comparisons(.x, data = attrition_treat_data), .id = "lhs", .progress = TRUE)

attrition_treat_comp_df %>%
  write_csv(file.path(script_options$output_path, "attrition_treat_comp_df.csv"))

# Joint tests (close, far, overall)
attrition_treat_joint_tests = map(
  attrition_treat_fits,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

list(
  joint_tests = attrition_treat_joint_tests,
  n_attrition = n_attrition
) %>%
  saveRDS(file.path(script_options$output_path, "attrition_treat_joint_tests.rds"))

# ---- Non-attrited HHs: are their characteristics balanced across treatment? ----
non_attrition_treat_data = attrition_endline_data %>%
  select(-where(is.list)) %>%
  filter(!not_in_know_table) %>%
  mutate(
    treat_dist = paste0("treat: ", assigned.treatment, ", dist: ", dist.pot.group) %>% factor()
  )

n_non_attrition = nrow(non_attrition_treat_data)

non_attrition_treat_fits = feols(
  data = non_attrition_treat_data,
  .[names(attrition_balance_vars)] ~ 0 + treat_dist + i(county, ref = "Busia"),
  cluster = ~cluster.id
)

non_attrition_treat_comp_df = non_attrition_treat_fits %>%
  map_dfr(~create_balance_comparisons(.x, data = non_attrition_treat_data), .id = "lhs", .progress = TRUE)

non_attrition_treat_comp_df %>%
  write_csv(file.path(script_options$output_path, "non_attrition_treat_comp_df.csv"))

non_attrition_treat_joint_tests = map(
  non_attrition_treat_fits,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

list(
  joint_tests = non_attrition_treat_joint_tests,
  n_non_attrition = n_non_attrition
) %>%
  saveRDS(file.path(script_options$output_path, "non_attrition_treat_joint_tests.rds"))

}

# ============================================================================
# Monitoring Attrition Balance Checks
# ============================================================================
# 991 wave 1 individuals were eligible for PoT monitoring (monitored == TRUE
# in the census) but didn't make it onto the monitoring roster (true.monitored
# == FALSE). These individuals are excluded from the structural model. We check
# whether this attrition is correlated with treatment.

if (run_all || script_options$monitored_attrition) {

library(fixest)

census_data = census_data %>%
  filter(!is.na(cluster.id)) 


# Construct survey_cto_dropped from census_data (has monitored & true.monitored)
# and join onto full_analysis_data for analysis variables


mon_attrition_data = census_data %>%
  filter(!is.na(cluster.id)) %>%
  filter(sms.treatment == "sms.control") %>%
  filter(monitored) %>%
  filter(have_phone == "No" | have_phone == "Don't know number") %>%
  mutate(
    survey_cto_dropped = monitored & !true.monitored
  ) %>%
  select(KEY.individ, survey_cto_dropped, monitored, true.monitored) %>%
  inner_join(
    full_analysis_data %>%
      select(KEY.individ, assigned.treatment, dist.pot.group, county, cluster.id,
             age.census, gender, have_phone, cluster.dist.to.pot, name_matched),
    by = "KEY.individ"
  ) %>%
  filter(!is.na(assigned.treatment)) %>%
  mutate(
    female = as.numeric(gender == "female"),
    have_phone_lgl = as.numeric(have_phone == "Yes")
  ) %>%
  left_join(
    n_indiv_df %>% transmute(cluster.id, n_per_cluster),
    by = "cluster.id"
  ) 




cat("\n--- Monitoring Attrition Summary ---\n")
cat("Total non-Phone SMS control individuals:", nrow(mon_attrition_data), "\n")
cat("Dropped from monitoring roster:", sum(mon_attrition_data$survey_cto_dropped), "\n")

# ---- Test 1: Does treatment predict survey_cto_dropped? ----
mon_attrition_pooled = mon_attrition_data %>%
  feols(
    survey_cto_dropped ~ i(assigned.treatment, ref = "control") | county,
    cluster = ~cluster.id
  )

mon_attrition_treat_dist = mon_attrition_data %>%
  feols(
    survey_cto_dropped ~ i(assigned.treatment, dist.pot.group, ref = "control") | county,
    cluster = ~cluster.id
  )

setFixest_dict(c(
  survey_cto_dropped = "Dropped from Monitoring",
  "assigned.treatment::ink" = "Ink",
  "assigned.treatment::calendar" = "Calendar",
  "assigned.treatment::bracelet" = "Bracelet"
))

mon_attrition_tex_postprocessing = function(tex) {
  tex %>%
    str_remove("\\\\begin\\{table\\}\\[htbp\\]") %>%
    str_remove("\\\\end\\{table\\}") %>%
    str_replace(., "Covariate", "\\\\midrule Covariate")
}

etable(mon_attrition_pooled, mon_attrition_treat_dist,
       headers = c("Pooled", "Treat $\\times$ Distance"),
       depvar = FALSE,
       fitstat = c("n", "r2"),
       se.below = TRUE,
       tex = TRUE,
       title = "Differential Monitoring Attrition by Treatment Assignment",
       label = "tab:monitoring-attrition-treatment",
       notes = "",
       file = "presentations/tables/monitoring-attrition-by-treatment.tex",
       replace = TRUE,
       postprocess.tex = mon_attrition_tex_postprocessing,
       digits = 3,
       digits.stats = 3,
       drop.section = "fixef",
       style.df = style.df(depvar.title = "", fixef.title = "", var.title = "", stats.title = "")
)

# ---- Test 2: Are census covariates balanced between kept and dropped? ----
mon_attrition_balance_vars = c(
  age.census = "Age",
  female = "Female",
  n_per_cluster = "Number of individuals per community",
  cluster.dist.to.pot = "Distance to PoT"
)

mon_attrition_balance_df = map_dfr(names(mon_attrition_balance_vars), function(v) {
  fit = feols(as.formula(paste0(v, " ~ survey_cto_dropped ")),
              data = mon_attrition_data, cluster = ~cluster.id)
  ct = coeftable(fit)

  group_means = mon_attrition_data %>%
    summarise(
      mean_present = mean(.data[[v]][!survey_cto_dropped], na.rm = TRUE),
      mean_missing = mean(.data[[v]][survey_cto_dropped], na.rm = TRUE),
      n_present = sum(!is.na(.data[[v]]) & !survey_cto_dropped),
      n_missing = sum(!is.na(.data[[v]]) & survey_cto_dropped)
    )

  tibble(
    variable = v,
    label = mon_attrition_balance_vars[v],
    mean_present = group_means$mean_present,
    mean_missing = group_means$mean_missing,
    diff = ct["survey_cto_droppedTRUE", "Estimate"],
    se = ct["survey_cto_droppedTRUE", "Std. Error"],
    pval = ct["survey_cto_droppedTRUE", "Pr(>|t|)"],
    n_present = group_means$n_present,
    n_missing = group_means$n_missing
  )
})

cat("\n--- Monitoring Attrition Covariate Balance ---\n")
mon_attrition_balance_df %>%
  mutate(stars = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**", pval < 0.1 ~ "*", TRUE ~ "")) %>%
  select(label, mean_present, mean_missing, diff, se, pval, stars) %>%
  print(n = Inf)

# ---- .tex output ----
mon_fmt = function(x, d = 3) formatC(x, format = "f", digits = d)
mon_stars_tex = function(p) case_when(p < 0.01 ~ "^{***}", p < 0.05 ~ "^{**}", p < 0.1 ~ "^{*}", TRUE ~ "")

mon_tex_rows = mon_attrition_balance_df %>%
  mutate(row = paste0(
    label, " & ",
    mon_fmt(mean_present), " & ",
    mon_fmt(mean_missing), " & ",
    "$", mon_fmt(diff), mon_stars_tex(pval), "$ & ",
    "$(", mon_fmt(se), ")$ & ",
    format(n_present + n_missing, big.mark = ","),
    " \\\\"
  )) %>%
  pull(row)

mon_n_present_total = max(mon_attrition_balance_df$n_present)
mon_n_missing_total = max(mon_attrition_balance_df$n_missing)

mon_tex_table = c(
  "\\centering",
  "\\begin{tabular}{l ccccc}",
  "\\hline\\hline",
  paste0(" & \\multicolumn{1}{c}{Monitored} & \\multicolumn{1}{c}{Dropped} & ",
         "\\multicolumn{1}{c}{Difference} & \\multicolumn{1}{c}{SE} & \\multicolumn{1}{c}{N} \\\\"),
  "\\hline",
  "\\addlinespace",
  "\\textit{Census Variables} \\\\",
  "\\addlinespace",
  mon_tex_rows,
  "\\hline\\hline",
  "\\end{tabular}"
)

writeLines(mon_tex_table, "presentations/tables/monitoring-attrition-covariate-balance.tex")

# ---- Test 3: Treatment balance among dropped individuals ----
mon_attrition_treat_data = mon_attrition_data %>%
  filter(survey_cto_dropped) %>%
  mutate(
    treat_dist = paste0("treat: ", assigned.treatment, ", dist: ", dist.pot.group) %>% factor()
  )

n_mon_attrition = nrow(mon_attrition_treat_data)

mon_attrition_treat_fits = feols(
  data = mon_attrition_treat_data,
  .[names(mon_attrition_balance_vars)] ~ 0 + treat_dist,
  cluster = ~cluster.id
)

# Pairwise comparisons
mon_attrition_treat_comp_df = mon_attrition_treat_fits %>%
  map_dfr(~create_balance_comparisons(.x, data = mon_attrition_treat_data), .id = "lhs", .progress = TRUE)

mon_attrition_treat_comp_df %>%
  write_csv(file.path(script_options$output_path, "monitoring_attrition_comp_df.csv"))

# Joint tests (close, far, overall)
mon_attrition_treat_joint_tests = map(
  mon_attrition_treat_fits,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

list(
  joint_tests = mon_attrition_treat_joint_tests,
  n_mon_attrition = n_mon_attrition
) %>%
  saveRDS(file.path(script_options$output_path, "monitoring_attrition_joint_tests.rds"))

# ---- Monitoring non-dropped: are their characteristics balanced across treatment? ----
mon_non_attrition_treat_data = mon_attrition_data %>%
  filter(!survey_cto_dropped) %>%
  mutate(
    treat_dist = paste0("treat: ", assigned.treatment, ", dist: ", dist.pot.group) %>% factor()
  )

n_mon_non_attrition = nrow(mon_non_attrition_treat_data)

mon_non_attrition_treat_fits = feols(
  data = mon_non_attrition_treat_data,
  .[names(mon_attrition_balance_vars)] ~ 0 + treat_dist,
  cluster = ~cluster.id
)

mon_non_attrition_treat_comp_df = mon_non_attrition_treat_fits %>%
  map_dfr(~create_balance_comparisons(.x, data = mon_non_attrition_treat_data), .id = "lhs", .progress = TRUE)

mon_non_attrition_treat_comp_df %>%
  write_csv(file.path(script_options$output_path, "mon_non_attrition_comp_df.csv"))

mon_non_attrition_treat_joint_tests = map(
  mon_non_attrition_treat_fits,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

list(
  joint_tests = mon_non_attrition_treat_joint_tests,
  n_mon_non_attrition = n_mon_non_attrition
) %>%
  saveRDS(file.path(script_options$output_path, "mon_non_attrition_joint_tests.rds"))

}


# ============================================================================
# SMS Balance Analysis
# ============================================================================
if (run_all || script_options$sms) {

stop()

census_data %>%
  count(cluster.id)
  analysis_data %>%
    count(cluster.id)

  sms_itt_sample_df = census_data %>%
    filter(have_phone == "Yes") %>%
    filter(sms.treatment == "social.info" | (sms.treatment == "sms.control" & sms.ctrl.sample.order == 1))  %>%
    mutate(
      sms_treatment = case_when(
        sms.treatment == "sms.control" ~ "smscontrol",
        sms.treatment %in% c("reminder.only", "social.info") ~ "smstreatment",
        TRUE ~ NA_character_
      )
    )  %>%
    left_join(
      full_analysis_data %>%
        select(cluster.id, cluster.dist.to.pot) %>% unique,
      by = "cluster.id"
    ) %>%
    filter(!is.na(cluster.id)) %>%
    mutate(cluster.id_fac = factor(cluster.id)) %>%
    left_join(
      analysis_data %>%
        mutate(
            treat_dist = paste0(
            "treat: ", 
            assigned.treatment,
            ", dist: ", dist.pot.group
            ) 
          )   %>%
        select(cluster.id, treat_dist) %>% unique(),
      by = c("cluster.id_fac" = "cluster.id")
    ) %>%
    mutate(gender = case_when(
      gender == 2 ~ "female",
      gender == 1 ~ "male",
      TRUE ~ NA_character_
    ))


  sms_sample_social_info_df = full_analysis_data %>%
    filter(sms.treatment == "social.info" | sms.treatment == "sms.control") %>%
    mutate(
      sms_treatment = case_when(
        sms.treatment == "sms.control" ~ "smscontrol",
        sms.treatment %in% c("reminder.only", "social.info") ~ "smstreatment",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(sms_treatment)) %>%
    filter(have_phone == "Yes") %>%
    # subset to SMS treated OR first SMS control - we shouldn't be adding extra 
    # SMS control HHs not in PAP
    filter(
        (sms_treatment == "smstreatment") |
        (sms_treatment == "smscontrol" & sms.ctrl.sample.order == 1)
    ) 

  sms_enrollment_balance_fit = sms_itt_sample_df %>%
    filter(sms_treatment == "smstreatment") %>%
    mutate(
      in_ana_df = KEY.individ %in% sms_sample_social_info_df$KEY.individ
    ) %>%
    feols(
      in_ana_df ~ 0 + treat_dist + i(county, ref = "Busia"),
      cluster = ~cluster.id
    )


  create_sms_bal_data = function(sms_df, n_indiv_df) {
    sms_df %>%
    select(
      sms_treatment, cluster.id, county,
      age.census, gender, have_phone, cluster.dist.to.pot) %>%
    left_join(
      n_indiv_df %>% transmute(cluster.id, n_per_cluster),
      by = "cluster.id"
    ) %>%
    transmute(
      sms_treatment = factor(sms_treatment, levels = c("smscontrol", "smstreatment")),
      cluster.id,
      county,
      age = age.census,
      female = case_when(
        gender == "female" ~ TRUE,
        gender == "male" ~ FALSE,
        TRUE ~ NA
      ),
      have_phone_lgl = case_when(
        have_phone == "Yes" ~ TRUE,
        have_phone == "No" ~ FALSE,
        TRUE ~ NA
      ),
      n_per_cluster,
      cluster.dist.to.pot = cluster.dist.to.pot / 1000
    )
  }


  sms_bal_data = create_sms_bal_data(sms_sample_social_info_df, n_indiv_df)
  sms_itt_df = create_sms_bal_data(sms_itt_sample_df, n_indiv_df)


  sms_bal_vars = c("age", "female", "n_per_cluster", "cluster.dist.to.pot")


  sms_itt_bal_fit = feols(
    data = sms_itt_df,
    .[sms_bal_vars] ~ 0 +   sms_treatment + i(county, ref = "Busia"),
    cluster = ~cluster.id,
    data.save = TRUE
  )
  sms_bal_fit = feols(
    data = sms_bal_data,
    .[sms_bal_vars] ~ 0+ sms_treatment + i(county, ref = "Busia"),
    cluster = ~cluster.id,
    data.save = TRUE
  )

  age_itt_df = sms_itt_df %>%
    filter(sms_treatment == "smscontrol") %>%
    select(age, sms_treatment, county) %>%
    arrange(age)

  age_bal_df = sms_bal_data %>%
    filter(sms_treatment == "smscontrol") %>%
    select(age, sms_treatment, county) %>%
    arrange(age)
get_control_mean_from_fit = function(fit,
                                     treatment_var = "sms_treatment",
                                     control_level = "smscontrol") {
  data = fit$data
  y_var = as.character(fit$fml_all$linear[[2]])
  outcome_y = data[[y_var]][data[[treatment_var]] == control_level]
  list(
    mu = mean(outcome_y, na.rm = TRUE),
    se = sd(outcome_y, na.rm = TRUE) / sqrt(sum(!is.na(outcome_y)))

  )
}
sms_comp_fn = function(fit) {
  tidy_fit = tidy(fit, conf.int = TRUE)

  control_row = tidy_fit %>%
    filter(str_detect(term, "smscontrol")) %>%
    transmute(
      lhs_treatment = "smscontrol",
      rhs_treatment = NA_character_,
      lhs_dist = "combined",
      rhs_dist = "combined",
      estimate, std.error, statistic, p.value,
      comp_type = "adjusted_control"
    )

  treat_comp = hypotheses(
    fit,
    "sms_treatmentsmstreatment - sms_treatmentsmscontrol = 0"
  ) %>%
    as_tibble() %>%
    transmute(
      lhs_treatment = "smstreatment",
      rhs_treatment = "smscontrol",
      lhs_dist = "combined",
      rhs_dist = "combined",
      estimate, std.error, statistic, p.value,
      comp_type = "treatment"
    )

  control_mean_row = tibble(
    term = "control_mean",
    lhs_treatment = "smscontrol",
    rhs_treatment = NA_character_,
    lhs_dist = "combined",
    rhs_dist = "combined",
    estimate = get_control_mean_from_fit(fit)$mu,
    std.error = get_control_mean_from_fit(fit)$se,
    statistic = NA_real_,
    p.value = NA_real_,
    comp_type = "control_mean"
  )

  bind_rows(control_row, treat_comp, control_mean_row)
}


  sms_comp_df = sms_bal_fit %>%
    map_dfr(
      sms_comp_fn,
      .id = "lhs"
    )
  sms_itt_comp_df = sms_itt_bal_fit %>%
    map_dfr(
      sms_comp_fn,
      .id = "lhs"
    )

  sms_comp_df %>%
    filter(lhs == "lhs: age") %>%
    filter(lhs_treatment == "smscontrol")

  sms_itt_comp_df %>%
    filter(lhs == "lhs: age") %>%
    filter(lhs_treatment == "smscontrol")



  N_sms_control = sum(sms_bal_data$sms_treatment == "smscontrol")
  N_sms_treat   = sum(sms_bal_data$sms_treatment == "smstreatment")
  N_sms_treat = 3022

  list(
    sms_comp_df    = sms_comp_df,
    N_sms_control  = N_sms_control,
    N_sms_treat    = N_sms_treat
  ) %>%
    saveRDS(file.path(script_options$output_path, "sms_balance.rds"))

  list(
    sms_itt_comp_df = sms_itt_comp_df,
    N_sms_control  = sum(sms_itt_df$sms_treatment == "smscontrol"),
    N_sms_treat    = sum(sms_itt_df$sms_treatment == "smstreatment")
  ) %>%
    saveRDS(file.path(script_options$output_path, "sms_itt_balance.rds"))

  sms_itt_comp_df %>%
    write_csv(file.path(script_options$output_path, "sms_itt_comp_df.csv"))

  sms_comp_df %>%
    write_csv(file.path(script_options$output_path, "sms_comp_df.csv"))

    # enrolment balance checks
  sms_enrollment_balance_results =
  perform_balance_joint_test(
   sms_enrollment_balance_fit,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )

comp_sms_enrollment = create_balance_comparisons(
  sms_enrollment_balance_fit,
  data = sms_itt_sample_df %>%
    filter(sms_treatment == "smstreatment") %>%
    mutate(
      in_ana_df = KEY.individ %in% sms_sample_social_info_df$KEY.individ
    )
)

  # save enrollment balance results
  list(
    enrollment_balance_comp = comp_sms_enrollment,
    enrollment_balance_joint_tests = sms_enrollment_balance_results
  ) %>%
    saveRDS(file.path(script_options$output_path, "sms_enrollment_balance.rds")
  )

}


# counts for mermaid diagram
# census_data = census_data %>%
#   filter(!is.na(cluster.id))


# census_data = census_data %>%
#   mutate(
#     non_na_endline_type = !is.na(endline.type),
#     sms_status = case_when(
#       sms.treatment == "sms.control" ~ "control",
#       sms.treatment %in% c("reminder.only", "social.info") ~ "treatment",
#       TRUE ~ NA_character_
#     )
#   )

# # SHould be SMS control 8,144, treat 3,022 
# census_data %>%
#   count(non_na_endline_type, sms_status)
# # 8,144

# census_data %>%
#   filter(sms.consent & true.monitored) %>%
#   count(non_na_endline_type, sms_status)
# #3,022

# # 11,166
# sub_df = census_data %>%
#   filter(
#     (sms.consent & true.monitored) | (non_census_data = census_data %>%
#   filter(!is.na(cluster.id))


# census_data = census_data %>%
#   mutate(
#     non_na_endline_type = !is.na(endline.type),
#     sms_status = case_when(
#       sms.treatment == "sms.control" ~ "control",
#       sms.treatment %in% c("reminder.only", "social.info") ~ "treatment",
#       TRUE ~ NA_character_
#     )
#   )

# # SHould be SMS control 8,144, treat 3,022 
# census_data %>%
#   count(non_na_endline_type, sms_status)
# # 8,144

# census_data %>%
#   filter(sms.consent & true.monitored) %>%
#   count(non_na_endline_type, sms_status)
# #3,022

# # 11,166
# sub_df = census_data %>%
#   filter(
#     (sms.consent & true.monitored) | (non_na_endline_type & sms.treatment == "sms.control") 
#   )
# # 2,990 control and 1,230 treat
# sub_df %>%
#   group_by(endline, sms_status) %>%
#   summarise(
#     n_key = n_distinct(KEY),
#     n_key_individ = n_distinct(KEY.individ)
#   )
# # 2,990 control and 1,228 actual treat


# # roster selected = 2990 + 1230 -> 4220

# colnames(sub_df)
# # Reached should be: 3830
# 2774 + 1056
# # nrow(all_endline_data) = 3,678
# nrow(all_endline_data)

# all_endline_data %>%
#   count(sms_status)

# all_endline_data %>%
#   filter(!is.na(cluster.id)) %>%
#   count(sms_status)
# na_endline_type & sms.treatment == "sms.control") 
#   )
# # 2,990 control and 1,230 treat
# sub_df %>%
#   group_by(endline, sms_status) %>%
#   summarise(
#     n_key = n_distinct(KEY),
#     n_key_individ = n_distinct(KEY.individ)
#   )
# # 2,990 control and 1,228 actual treat


# # roster selected = 2990 + 1230 -> 4220

# colnames(sub_df)
# # Reached should be: 3830
# 2774 + 1056
# # nrow(all_endline_data) = 3,678
# nrow(all_endline_data)

# all_endline_data %>%
#   count(sms_status)

# all_endline_data %>%
#   filter(!is.na(cluster.id)) %>%
#   count(sms_status)
