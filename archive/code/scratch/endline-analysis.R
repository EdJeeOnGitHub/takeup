library(tidyverse)
library(broom)
library(data.table)
library(kableExtra)
library(knitr)
library(fixest)



###### Structure ######
## Loading scripts and cleaned data
## Analysis of endline data in order:
# 1. Change in externality knowledge baseline vs endline
# 2. Endline externality knowledge regressions
# 3. Generating know/don't know beliefs
# 4. Endline beliefs regressions
# 5. Predicted deworming at endline
# 6. Incentive implementation checks
# 7. Gift preferences 
# 8. Travel time analysis
###### End Structure ######

##### 0. Setup -----------------------------------------------------------------
# Hyper parameters for saving output/config options
params = lst(
    table_output_path = "presentations/rf-tables/main-specs",
    show_probs = FALSE,
    width = 0.95,
    cache = FALSE,
    fit = FALSE,
    stat = "std.error" # "ci", "p", "std.error"
)
# Karim's helper functions
source(file.path("R/common/analysis.R"))
# Load reduced form data
source(file.path("scratch", "reduced-form-setup.R"))

# LASSO Covariates
# From running:
# pdslasso dewormed_num dpf ($cov_vars i.county_fac mu_d), cluster(clusteridx) pnotpen(i.county_fac)
# where mu_d is the expected distance to the cluster
# get the same result using actual distance.
l_cov_vars = c(
  "female",
  "age.census"
)

###### 1. Change in externality knowledge baseline vs endline ------------------
#### Change in Externality Knowledge (Table 13)
# Checking endline vs baseline externality knowledge
balance_data = read_rds(
  file.path(
    "temp-data",
    "saved_balance_data.rds"
  )
)

balance_analysis_data = balance_data$analysis_data
baseline_worm_data = balance_data$baseline_worm_data
endline_vars = balance_data$endline_vars
endline_and_baseline_worm_data = balance_data$endline_and_baseline_worm_data
pretreat_data = balance_data$pretreat_data
clean_census_data = balance_data$clean_census_data

# Function to create endline vs baseline treatment effect tables
create_te_table = function(data, var) {
  data %>%
  group_by(
    treat_dist,
    type
  ) %>%
  summarise(
    mean = mean({{var}}, na.rm = TRUE)
  ) %>%
  mutate(
    treat = str_extract(treat_dist, "(?<=treat\\: )\\w+"),
    dist = str_extract(treat_dist, "(?<=dist\\: )\\w+"),
  ) %>%
  pivot_wider(names_from = type, values_from = mean)  %>%
  group_by(dist) %>%
  mutate(
    baseline_te = baseline - baseline[treat == "control"],
    endline_te = endline - endline[treat == "control"]
  ) %>%
  select(-treat_dist)
}


externality_data = endline_data %>%
    mutate(
      fully_aware_externalities = case_when(
        neighbours_worms_affect == "yes" & worms_affect == "yes" ~ TRUE, 
        # Ed: 2025-08-08 NA in these two variables is actually "don't know" due to 
        # a coding error in `R/common/analysis.R:129` in SurveyCTO these two 
        # variables use different binary encoding for yes/no and the original 
        # code corrects this but doesn't correct "don't know" correctly
        is.na(neighbours_worms_affect) | is.na(worms_affect) ~ FALSE,
        TRUE ~ FALSE
      ),
      know_worms_infectious = spread_worms == "yes",
      externality_omnibus = fully_aware_externalities | know_worms_infectious
    ) %>%
    select(KEY.individ, cluster.id, externality_omnibus) 

# Just the baseline data, with cluster level covariates joined on
# Not used in paper, but for checking
solo_baseline_data = baseline_data  %>%
    transmute(
      fully_aware_externalities = case_when(
        neighbours_worms_affect == "yes" & worms_affect == "yes" ~ TRUE, 
        is.na(neighbours_worms_affect) | is.na(worms_affect) ~ NA,
        TRUE ~ FALSE
      ),
      know_worms_infectious = spread_worms == "yes",
      externality_omnibus = fully_aware_externalities | know_worms_infectious,
      cluster.id
    ) %>%
    left_join(
      cov_analysis_data %>%
        select(
          standard_cluster.dist.to.pot,
          assigned_treatment,
          assigned_dist_group,
          mu_d,
          cluster_id
        ) %>% unique(),
        by = c("cluster.id" = "cluster_id")
    )

# Just the endline data, with cluster level covariates joined on
# Not used in paper, but for checking
solo_endline_data = endline_data %>%
    transmute(
      fully_aware_externalities = case_when(
        neighbours_worms_affect == "yes" & worms_affect == "yes" ~ TRUE, 
        is.na(neighbours_worms_affect) | is.na(worms_affect) ~ NA,
        TRUE ~ FALSE
      ),
      know_worms_infectious = spread_worms == "yes",
      externality_omnibus = fully_aware_externalities | know_worms_infectious,
      cluster.id
    ) %>%
    left_join(
      cov_analysis_data %>%
        select(
          standard_cluster.dist.to.pot,
          assigned_treatment,
          assigned_dist_group,
          mu_d,
          cluster.id.x
        ) %>% unique(),
        by = c("cluster.id" = "cluster.id.x")
    )

# Comparison of endline vs baseline knowledge by treatment/dist group
# Not used in paper, but for checking
solo_comp_df = inner_join(
  solo_endline_data %>%
    group_by(assigned_treatment, assigned_dist_group) %>%
    summarise(
      mean_know_endline = mean(know_worms_infectious, na.rm = TRUE)
    ),
  solo_baseline_data %>%
    group_by(assigned_treatment, assigned_dist_group) %>%
    summarise(
      mean_know_baseline = mean(know_worms_infectious, na.rm = TRUE)
    ),
    by = c("assigned_treatment", "assigned_dist_group")
)

###### 2. Endline externality knowledge regressions ----------------------------
# Externality knowledge data with cluster level variables + indiv covariates 
externality_knowledge_df = cov_analysis_data %>%
  select(
    cluster_id, 
    cluster.id,
    cluster.id.x,
    assigned_treatment,
    assigned_dist_group,
    mu_d,
    standard_cluster.dist.to.pot,
    county,
    all_of(l_cov_vars),
    KEY.individ
    )  %>%
    inner_join(
      externality_data %>% select(-cluster.id),
      by = "KEY.individ"
    ) %>%
    mutate(cluster_id_rank = dense_rank(cluster.id)) 



# Wider sample externality knowledge data, not main analysis sample
full_externality_knowledge_df = full_analysis_data  %>%
  mutate(
    female = gender == "female",
    cluster_id = dense_rank(cluster.id)
    ) %>%
  select(
    cluster_id,
    cluster.id,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group,
    county,
    all_of(l_cov_vars),
    KEY.individ,
    cluster_id_rank,
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

externality_knowledge_df %>% nrow()
full_externality_knowledge_df %>% nrow()

# Externality knowledge regression function
externality_knowledge_regression = function(data, weights) {
  feols(
    externality_omnibus ~ 0 + assigned_treatment*assigned_dist_group + .[l_cov_vars] + mu_d   | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}


# Wrapper function runs regression, calculates SEs, and saves output
externality_knowledge_output = wrapper_function(
  data = externality_knowledge_df,
  regression_spec = externality_knowledge_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/externality-knowledge-tidy-tes.csv",
  table_name = "rf_externality_knowledge_tbl",
  table_options = list(
    dependent_var = "Dependent variable: Externality Knowledge"
  )
)
# robustness to including non-monitored individuals in the 
# externality knowledge regression
full_externality_knowledge_output = wrapper_function(
  data = full_externality_knowledge_df %>%
    filter(sms.treatment.2 == "sms.control"),
  regression_spec = externality_knowledge_regression,
  tidy_summ_path = "temp-data/tidy-rf-tes/full-externality-knowledge-tidy-tes.csv",
  table_name = "rf_full_externality_knowledge_tbl",
  table_options = list(
    dependent_var = "Dependent variable: Externality Knowledge"
  )
)



###### 3. Generating know/don't know beliefs -----------------------------------

# Summary table of beliefs by relationship type
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

# Beliefs data for main analysis sample
belief_ana_df = analysis_data %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  nest_join(
    endline_know_table_data %>% 
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

# Joining on cluster level distance variables
all_data = analysis.data %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


# Disaggregated belief data for main analysis sample + wider sample
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
    endline_know_table_data %>% 
      filter(fct_match(know.table.type, "table.A")) %>%
      mutate(ed_in_endline_flag = TRUE),
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
        ed_in_endline_flag = any(.x$ed_in_endline_flag, na.rm = TRUE)
      )
    }
  )) %>%
    # filter(obs_know_person > 0)  %>%
    select(
      KEY.individ, 
      contains("know"), 
      assigned.treatment, 
      sms.treatment,
      ed_in_endline_flag,
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
           sms.treatment,
           obs_know_person,
           ed_in_endline_flag,
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


# Calculating how many people recognize at least one person
# per sub-sample in our full belief data
disagg_belief_all_df %>%
  filter(ed_in_endline_flag == TRUE) %>%
  group_by(
    not_mon = is.na(dewormed), 
    sms.treatment,
    obs_at_least_1 = obs_know_person > 0
    ) %>%
  summarise(
    n = n_distinct(KEY.individ)
  ) %>%
  filter(
    (not_mon == FALSE & sms.treatment == "reminder.only") |
    (not_mon == FALSE & sms.treatment == "sms.control") |
    (not_mon == FALSE & sms.treatment == "social.info") |
    (not_mon == TRUE & sms.treatment == "sms.control")
  ) %>%
  pivot_wider(
    names_from = obs_at_least_1,
    names_prefix = "obs_at_least_1_",
    values_from = n
  ) %>%
  rename(
    recognize_0 = obs_at_least_1_FALSE,
    recognize_someone = obs_at_least_1_TRUE
  ) 

 disagg_belief_all_df = disagg_belief_all_df %>%
  filter(obs_know_person > 0) 


# Disaggregated belief data for main analysis sample
disagg_base_belief_data = cov_analysis_data %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  nest_join(
    endline_know_table_data %>% 
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


# Taking disagg belief data for main analysis sample and converting to a
# know/don't know df
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

# Taking disagg belief data for wider sample and converting to a know/don't know
# df
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

# First and second order belief dfs
know_1_df = know_df  %>%
  filter(belief_type == "1ord") 
know_2_df = know_df  %>%
  filter(belief_type == "2ord")

# first order beliefs for full sample
know_1_all_df = know_all_df  %>%
  filter(belief_type == "1ord") %>%
  mutate(cluster_id = as.numeric(cluster.id)) 

###### 4. Endline beliefs regressions ------------------------------------------
# Belief regression functions
## Discrete
discrete_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + assigned_dist_group + i(assigned_treatment, assigned_dist_group, "control") + .[l_cov_vars] +  mu_d | county,
    data = data,
    weights = weights
  )
}
## Continuous
cts_f_know = function(data, weights) {
  feols(
    prop_knows ~ assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") + .[l_cov_vars] + mu_d | county,
    data = data,
    weights = weights
  )
}
## Household Distance
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

# No meaningful difference
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


###### 5. Predicted deworming at endline ---------------------------------------
####  Endline Predicted Deworming Takeup
endline_data_full = endline.data %>%
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
    analysis.data %>%
      select(KEY.individ, gender, age.census, cluster.id) %>%
      mutate(female = gender == "female")  %>%
      left_join(
        cov_analysis_data %>%
          select(cluster.id.x, mu_d) %>%
          unique(),
          by = c("cluster.id" = "cluster.id.x")
      ) %>%
      select(KEY.individ, female, age.census, mu_d),
      by = "KEY.individ"
  )

endline_data = endline_data_full %>%
  filter(true.monitored == TRUE & sms.treatment == "sms.control")


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
endline_data_full %>%
    filter(sms.treatment == "sms.control" & true.monitored == TRUE)

# monitored + non monitored sample - but still SMS control
pred_dworm_full_output = wrapper_function(
  data = endline_data_full %>%
    filter(sms.treatment == "sms.control"),
  regression_spec = pred_dworm_fit,
  tidy_summ_path = "temp-data/tidy-rf-tes/predicted-endline-deworm-takeup-full-sample-tidy-tes.csv",
  table_name = "predicted_endline_deworm_takeup_full_spec_tbl",
  table_options = list(
    caption = "Average Treatment Effects: Reduced Form", 
    dependent_var = "Dependent variable: Predicted Take-up", 
    type = "APE", 
    stars = TRUE,
    drop_H0s = TRUE
    )
)

pred_dworm_full_output$tidy_summary %>%
  print(n = 100)


pred_dworm_full_output$tidy_summary  %>%
  filter(assigned_treatment == "bracelet")
pred_dworm_output$tidy_summary  %>%
  filter(assigned_treatment == "bracelet")


###### 6. Incentive implementation checks --------------------------------------
# Incentive Implementation 
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

###### 7. Gift preferences  ----------------------------------------------------
# Preference for Gift Fit Not Dewormed 
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
    group_by(assigned.treatment, dist.pot.group, dewormed) %>% 
    mutate(arm.size = n()) %>% 
    group_by(gift_choice, add = TRUE) %>%
    ungroup() %>%
    select(KEY.individ, cluster.id, gift_choice, 
      assigned.treatment, dist.pot.group, county, gender,
      age.census,
      cluster_id_rank
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
      )

  analysis.data  %>%
  # 38,019
    filter(!is.na(gift_choice)) %>%
    # 3,676
    filter(monitored) %>%
    # 3,329
    filter(monitor.consent) %>%
    # 3,329 %>%
    filter(dewormed == FALSE)
    # 1,808

  analysis.data  %>%
  # 38,019
    filter(!is.na(gift_choice)) %>%
    # 3,676
    filter(monitored) %>%
    # 3,329
    filter(monitor.consent)  %>%
    # 3,329
    filter(!hh.baseline.sample.pool)  %>%
    # 2,820
    filter(dewormed == FALSE)
    # 1,566


  analysis.data  %>%
  # 38,019
    filter(!is.na(gift_choice)) %>%
    # 3,676
    filter(monitored) %>%
    # 3,329
    filter(monitor.consent)  %>%
    # 3,329
    filter(!hh.baseline.sample.pool)  %>%
    # 2,820
    filter(sms.treatment.2 == "sms.control") %>%
    # 1,940
    filter(dewormed == FALSE)
    # 1,174


analysis_data %>%
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
    )  


summ_gift_df = analysis.data %>%
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
    select(KEY.individ, cluster.id, gift_choice, assigned.treatment,
           dist.pot.group, county, standard_cluster.dist.to.pot, cluster_id_rank) %>%
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





###### 8. Travel time analysis -------------------------------------------------
# Travel Time 
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




