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
  "lhs: pct_avail_message", "CHV share where to get dewormed?", "implementation"
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

analysis_data %>%
  select(cluster.id)

census_data %>%
  select(cluster.id)

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

# Adding school (PoT) data to analysis df
analysis_data = analysis_data %>%
    group_by(cluster.id) %>%
    mutate(row_id = 1:n()) %>%
    left_join(
      n_indiv_df %>% mutate(row_id = 1, cluster.id = factor(cluster.id)),
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

baseline_data


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

# Removing list columns as they mess up marginaleffects
endline_worm_data = endline_worm_data %>%
  select(where(~!is.list(.)))
endline_implementation_data = endline_implementation_data %>%
  select(where(~!is.list(.)))
baseline_worm_data = baseline_worm_data %>%
  select(where(~!is.list(.)))
pretreat_data = pretreat_data %>%
  select(where(~!is.list(.)))
clean_sens_imp_df = clean_sens_imp_df %>%
  select(where(~!is.list(.)))
praise_df = praise_df %>%
  select(where(~!is.list(.)))
stigma_df = stigma_df %>%
  select(where(~!is.list(.)))
pretreat_data = pretreat_data %>%
  select(where(~!is.list(.)))
clean_census_data = clean_census_data %>%
  select(where(~!is.list(.)))
analysis_data = analysis_data %>%
  select(where(~!is.list(.)))


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
      left_join(n_indiv_df %>% transmute(cluster.id = factor(cluster.id), n_per_cluster), by = "cluster.id") %>%
      select(any_of(takeup_vars), treat_dist, county, cluster.id) %>%
      # convert to Km for this table
      mutate(cluster.dist.to.pot = cluster.dist.to.pot / 1000 ), 
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


#### Save data so we can use in PDS estimates ----------------------------------
# takeup_vars
clust_n_df = n_indiv_df %>%
  select(cluster.id, n_per_cluster)
# pretreatment vars - pretreat_vars
clust_pretreat_df = pretreat_data %>%
  select(cluster.id, all_of(pretreat_vars)) %>%
  group_by(cluster.id) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
# baseline knowledge  - worm_vars
clust_baseline_worm_df = baseline_worm_data %>%
  select(cluster.id, all_of(worm_vars)) %>%
  group_by(cluster.id) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
# social image concerns  - praise_vars & stigma_vars
clust_praise_df = praise_df %>%
  select(cluster.id, all_of(praise_vars)) %>%
  group_by(cluster.id) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
clust_stigma_df = stigma_df %>%
  select(cluster.id, all_of(stigma_vars)) %>%
  group_by(cluster.id) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

cov_vars = c(
  "n_per_cluster",
  pretreat_vars,
  worm_vars,
  praise_vars,
  stigma_vars,
  takeup_vars
)

covariates_we_want = 
lhs_translation_df %>%
  mutate(
    lhs = str_remove(lhs, "lhs: ")
  ) %>%
  pull(lhs)

cov_vars = intersect(cov_vars, covariates_we_want)

covariate_df = analysis_data %>%
  mutate(cluster.id = as.numeric(levels(cluster.id)[cluster.id])) %>%
  select(cluster.id, any_of(takeup_vars), KEY.individ) %>%
  left_join(
    clust_n_df,
    by = "cluster.id"
  ) %>%
  left_join(
    clust_pretreat_df,
    by = "cluster.id"
  ) %>%
  left_join(
    clust_baseline_worm_df,
    by = "cluster.id"
  ) %>%
  left_join(
    clust_praise_df,
    by = "cluster.id"
  ) %>%
  left_join(
    clust_stigma_df,
    by = "cluster.id"
  ) %>%
  select(
    cluster.id, KEY.individ, all_of(cov_vars)
  )


# analysis_cov_df = analysis_data %>%
#   select(-all_of(takeup_vars)) %>%
#   left_join(
#     covariate_df,
#     by = "KEY.individ"
#   ) 


# analysis_cov_df %>%
#   write_csv(
#     "temp-data/analysis-cluster-covariate-data.csv"
#   )




create_balance_comparisons = function(fit) {
  comp_df = avg_comparisons(
    fit,
    variables = list("treat_dist" = "all")
  ) %>%
  as_tibble()


  comp_df = comp_df %>%
    mutate(
      lhs_treatment = str_extract(contrast, "treat: \\w+") %>% str_remove("treat: "),
      # remove first treat word and search after it for next treat word
      sub_str = str_extract(contrast, "(?<=treat).*"),
      rhs_treatment = str_extract(sub_str, "treat: \\w+") %>% str_remove("treat: "),

      lhs_dist = str_extract(contrast, "(?<=dist: )\\w+"),
      sub_str_dist = str_extract(contrast, "(?<=dist: ).*"),
      rhs_dist = str_extract(sub_str_dist, "(?<=dist: )\\w+")
    )  %>%
    select(-sub_str, -sub_str_dist)



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
    sample_mean_df = tibble(
      lhs_treatment = "control",
      rhs_treatment = NA,
      lhs_dist = "combined",
      rhs_dist = "combined",
      estimate = fitstat(fit, type = "my", verbose = FALSE)$my
    )

  rearranged_comp_df = rearranged_comp_df %>%
    bind_rows(
      control_mean_df,
      sample_mean_df
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
    .id = "lhs",
    .progress = TRUE
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

gen_close_p_val = function(x){
  map(
      c(split(worm_joint_hyp_matrix[1:4, ], 1:4), list(worm_joint_hyp_matrix[1:4, ])),
      ~perform_worm_change_test(endline_and_baseline_worm_data, x, .x)
  )
}

gen_far_p_val = function(x) {
  map(
    c(split(worm_joint_hyp_matrix[5:8, ], 1:4), list(worm_joint_hyp_matrix[5:8, ])),
    ~perform_worm_change_test(endline_and_baseline_worm_data,x,  .x)
)
}


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

takeup_dist_fit = feols(
  data = analysis_data %>%
    left_join(n_indiv_df %>% transmute(cluster.id = factor(cluster.id), n_per_cluster), by = "cluster.id") %>%
    select(any_of(takeup_vars), treat_dist, county, cluster.id), 
  .[takeup_vars] ~ 0 + cluster.dist.to.pot + i(county, ref = "Busia"),
  ~cluster.id
)



resample_cluster_dists = function(data, seed) {
  set.seed(seed)
  clust_dist_df = data %>%
    select(cluster.id, cluster.dist.to.pot) %>%
    unique()
  perm_clust_dist_df = clust_dist_df %>%
    mutate(
      perm_dist = sample(cluster.dist.to.pot, size = n())
    ) %>%
    select(-cluster.dist.to.pot)

    data = data %>%
      left_join(
        perm_clust_dist_df,
        by = "cluster.id"
      )
    return(data)
}


ri_fun = function(draw) {
  set.seed(draw)
  perm_baseline_worm_data = baseline_worm_data %>%
    resample_cluster_dists(draw)
  perm_census_data = clean_census_data %>%
    resample_cluster_dists(draw)

  perm_pretreat_data = pretreat_data %>%
    resample_cluster_dists(draw)

  perm_takeup_data = analysis_data %>%
    select(any_of(takeup_vars), cluster.dist.to.pot, county, cluster.id) %>%
    resample_cluster_dists(draw)

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

  perm_takeup_dist_fit = feols(
    data = perm_takeup_data, 
    .[takeup_vars] ~ 0 + perm_dist + i(county, ref = "Busia"),
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

analysis_data = analysis_data %>%
  left_join(
    n_indiv_df %>% transmute(cluster.id = factor(cluster.id), n_per_cluster), by = "cluster.id"
  )

if (script_options$fit_ri) {
  plan(multisession, workers = 12)
  perm_fit_df = future_map_dfr(
    1:500, 
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
    p_val = paste0("p = ", round(mean(statistic > realised_statistic), 3)),
    realised_statistic = unique(realised_statistic),
    x = quantile(statistic, 0.95, na.rm = TRUE),
    clean_name = unique(clean_name)
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
  ri_p_val_df,
  realised_fit_df
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

