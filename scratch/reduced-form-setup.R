
if (interactive()) {
library(tidyverse)
library(broom)
library(data.table)
library(kableExtra)
library(knitr)
library(fixest)
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
}

# Useful variables/hyperparameters
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
source(file.path("scratch", "reduced-form-functions.R"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"
## Loading Data - these are all intermediate cleaned data files
rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))
load(file.path("data", "takeup_village_pot_dist.RData"))
load(file.path("data", "analysis.RData"))


endline_data = ed_edits_prepare.endline.data(all.endline.data, census.data, cluster.strat.data)
# Remove endline.data - this is subset to a subsample we don't actually want
# instead, use all.endline.data
rm(endline.data)

baseline.data = read_rds("temp-data/reclean_baseline_data.rds") # Not sampling data!

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)



# Only take sms.control HHs
nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


# Only take monitored, no sms HHs
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

# main analysis sample
analysis_data <- monitored_nosms_data


# Save cluster distance sd in case we need to unstandardize - two variables 
# as there are two different naming conventions in use
sd_of_dist = sd(analysis_data$cluster.dist.to.pot)
dist_sd = sd(analysis_data$cluster.dist.to.pot)


## Load Census Data
# Setting up environment to avoid variable name clashes - this .RData file has 
# multiple datasets inside it with common df names
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

# Cleaning covariate data functions
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

clean_takeup_variables = function(data) {
  data %>%
    mutate(
      female = gender == "female",
      have_phone_lgl = phone_owner
    )
}


# Apply cleaning functions



## Cleaning baseline worm covariates
baseline_worm = baseline.data %>%
  clean_worm_covariates()
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

pretreat_data = clean_pretreat_covariates(baseline.data, endline_data) %>%
  left_join(cluster_treat_df, by = "cluster.id") %>%
  filter(!is.na(treat_dist))


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

# not sure what this does - inherited from KN
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

# Same with this
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


# Adding on Borusyak and Hull variables
# load expected distances
cell_expected_dist_df = read_csv(here::here("data", "cell_expected_dist_df.csv")) %>%
  rename(cell_expected_dist = dist, assigned_treatment = random_treat, assigned_dist_group = random_dist_group) 
cluster_expected_dist_df = read_csv(here::here("data", "cluster_expected_dist.csv")) %>%
  rename(clust_expected_dist = dist) %>%
  mutate(cluster.id = factor(cluster.id))


# Cleaning up analysis data coding of factors + adding BH variables
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

## Generating another version of analysis data with covariate data which we used to use 
# stata briefly
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



all_data = analysis.data %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

#### Create endline knowledge tables
raw_endline_know_table_data = bind_rows(table.A = read_csv(file.path("data", "raw-data", "Endline Survey-survey-sec_D-tableA.csv")) %>% 
                                       transmute(num.recognized = recogniseA,
                                                 dewormed = dewormedA,
                                                 second.order = order_2ndA,
                                                 second.order.reason = order_2nd_reason,
                                                 relationship = relationshipA,
                                                 relationship.other = relationshipA_other,
                                                 times.seen = times_seenA,
                                                 visited = visitedA, 
                                                 know.other.index = instanceA,
                                                 PARENT_KEY, KEY),
                                     table.B = read_csv(file.path("data", "raw-data", "Endline Survey-survey-sec_D-tableB.csv")) %>% 
                                       transmute(num.recognized = recogniseB,
                                                 dewormed = dewormedB,
                                                 dewormed.know.only = dewormedBB,
                                                 know.other.index = instanceB + 0:9, 
                                                 know.other.index.2 = know.other.index + 1, 
                                                 PARENT_KEY, KEY),
                                     .id = "know.table.type") %>% 
  mutate(recognized = num.recognized > 0)

raw_endline_know_table_data %>%
  summarise(
    n = n(),
    n_indiv = n_distinct(PARENT_KEY)
  )

survey_know_list = file.path("instruments", "SurveyCTO Forms", "Endline Survey", "Deployed Form Version", 
                              c("knowledge_list.csv", "kak_knowledge_list.csv")) %>% 
  map_dfr(read_csv, .id = "wave") %>% 
  mutate(wave = as.integer(wave))  %>%
  filter(wave == 1 | !is.na(survey.type))



endline_know_table_data = raw_endline_know_table_data %>%
  # Join on actual individuals from the endline
  left_join(select(endline_data, KEY.individ, wave, cluster.id, assigned.treatment, sms.treatment, dist.pot.group, KEY), 
             by = c("PARENT_KEY" = "KEY")) %>% 
  # Use the knowledge list so we know who knows each other
  left_join(select(survey_know_list, person_key, know.other.index, KEY.individ.other),
             by = c("KEY.individ" = "person_key", "know.other.index")) %>% 
  left_join(select(survey_know_list, person_key, know.other.index, KEY.individ.other), # Get the second person (table B)
            by = c("KEY.individ" = "person_key", "know.other.index.2" = "know.other.index"),
            suffix = c(".1", ".2")) %>% 
  mutate(across(c(recognized, dewormed, dewormed.know.only), factor, levels = c(0:1, 98), labels = c("no", "yes", "don't know"))) %>% 
  mutate(across(second.order, factor, levels = c(1:2, 97:98), labels = c("yes", "no", "prefer not say", "don't know"))) %>% 
  mutate(across(relationship, factor, levels = c(1:5, 99), labels = c("hh member", "extended family", "friend", "neighbor", "church", "other"))) %>%
  mutate(relationship = if_else(fct_match(relationship, "other") & fct_match(str_to_lower(relationship.other), c("village member", "village mate", "same village", "village elder", "village mates", "villager")),
                                "village member", as.character(relationship)) %>% factor,
         dewormed = if_else(fct_match(know.table.type, "table.B") & is.na(dewormed), dewormed.know.only, dewormed)) %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, respondent.dewormed.any = dewormed.any),
            by = "KEY.individ") %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, actual.other.dewormed.any.1 = dewormed.any),
            by = c("KEY.individ.other.1" = "KEY.individ")) %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, actual.other.dewormed.any.2 = dewormed.any),
            by = c("KEY.individ.other.2" = "KEY.individ")) %>% 
  mutate(actual.other.dewormed.any.either = actual.other.dewormed.any.1 | (fct_match(know.table.type, "table.B") & actual.other.dewormed.any.2)) 

# Missing some people in here - let's find from where
endline_know_table_data %>%
  filter(is.na(KEY.individ)) %>%
  summarise(
    n_parents = n_distinct(PARENT_KEY)
  )

endline_know_table_data %>%
  skimr::skim()



# HHs in endline, not in know table
in_endline_not_know_table = endline_data %>% 
  select(KEY.individ, KEY) %>%
  anti_join(raw_endline_know_table_data, by = c("KEY" = "PARENT_KEY")) %>%
  pull(KEY) %>%
  unique()

in_know_table_not_endline = raw_endline_know_table_data %>% 
  select(PARENT_KEY) %>%
  anti_join(endline_data, by = c("PARENT_KEY" = "KEY")) %>%
  pull(PARENT_KEY) %>%
  unique()

raw_endline_know_table_data %>%
  select(PARENT_KEY) %>%
  anti_join(all.endline.data %>% select(KEY), by = c("PARENT_KEY" = "KEY")) %>%
  pull(PARENT_KEY) %>%
  unique()  %>%
  length()

n_in_know_table = raw_endline_know_table_data %>%
  summarise(n = n_distinct(PARENT_KEY)) %>%
  pull(n)
n_in_endline = endline_data %>%
  summarise(n = n_distinct(KEY)) %>%
  pull(n)

n_in_know_table
n_in_endline
in_endline_not_know_table %>% length()
in_know_table_not_endline %>% length()


endline_know_table_data %>%
  summarise(
    n_parent_keys = n_distinct(PARENT_KEY),
    n_indiv_keys = n_distinct(KEY.individ),
    n_indiv_keys_know_treat = n_distinct(KEY.individ[!is.na(assigned.treatment)]),
    n_indiv_keys_know_dewormed = n_distinct(KEY.individ[!is.na(dewormed)])
  )

endline_data = endline_data %>%
  mutate(
    in_know_table = KEY %in% raw_endline_know_table_data$PARENT_KEY
  )

endline_data %>%
  group_by(in_know_table, dist.pot.group) %>%
  summarise(
    n = n()
  ) %>%
  mutate(pct = n / sum(n))

endline_data %>%
  group_by(in_know_table, assigned.treatment) %>%
  summarise(
    n = n()
  ) %>%
  mutate(pct = n / sum(n))

endline_data %>%
  group_by(in_know_table) %>%
  summarise(
    n = n(),
    n_monitored_consent = sum(monitor.consent == "yes", na.rm = TRUE),
    n_monitored = sum(monitored, na.rm = TRUE),
    n_reconsent = sum(reconsent),
  )


summ_endline_know_table = endline_know_table_data %>%
  filter(!is.na(KEY.individ)) %>%
  group_by(KEY.individ, know.table.type) %>%
  summarise(
    obs_know_person = sum(num.recognized),
    obs_know_person_prop = mean(num.recognized),
    knows_other_dewormed = sum(fct_match(dewormed, c("yes", "no")), na.rm = TRUE),
    knows_other_dewormed_yes = sum(fct_match(dewormed, "yes"), na.rm = TRUE),
    knows_other_dewormed_no = sum(fct_match(dewormed, "no"), na.rm = TRUE),
    thinks_other_knows = sum(fct_match(second.order, c("yes", "no")), na.rm = TRUE),
    thinks_other_knows_yes = sum(fct_match(second.order, "yes"), na.rm = TRUE),
    thinks_other_knows_no = sum(fct_match(second.order, "no"), na.rm = TRUE) 
  ) %>%
  ungroup()
summ_know_A_df = summ_endline_know_table %>%
  filter(know.table.type == "table.A") 

summ_know_B_df = summ_endline_know_table %>%
  filter(know.table.type == "table.B")

2686 + 1043

summ_endline_know_table %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )

endline_know_table_data %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )

summ_endline_know_table %>%
  left_join(
    analysis_data %>%
      select(KEY.individ, sms.treatment)  %>% 
      mutate(
        sms_status = case_when(
          sms.treatment == "sms.control" ~ "control", 
          sms.treatment %in% c("reminder.only", "social.info") ~ "treat",
          TRUE ~ NA_character_
        )
      ),
      by = "KEY.individ"
  ) %>%
  ungroup() %>%
  count(sms_status)


# Saving datasets to csv for debugging purposes - these are the datasets we use
# for the main analysis, but it's useful to have them in csv format to inspect
# and share with others more easily.

all.endline.data %>%
  write_csv(
    "temp-data/debug-all-endline-data.csv"
  )

analysis.data %>%
  write_csv(
    "temp-data/debug-analysis-data.csv"
  )

analysis_data %>%
  write_csv(
    "temp-data/analysis-data.csv"
  )