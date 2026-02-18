library(tidyverse)




# Village data
load(file.path("data", "takeup_village_pot_dist.RData"))
# Analysis data
load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

analysis.data = analysis.data %>%
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") 

# Only take sms.control HHs
nosms_data = analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


# Only take monitored, no sms HHs
monitored_nosms_data = analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


nosms_data %>%
  write_csv("data/clean-data/nosms-takeup-data.csv")

monitored_nosms_data %>%
  write_csv("data/clean-data/monitored-nosms-takeup-data.csv")


nosms_data %>%
  write_rds("data/clean-data/nosms-takeup-data.rds")

monitored_nosms_data %>%
  write_rds("data/clean-data/monitored-nosms-takeup-data.rds")

analysis.data %>%
  write_rds("data/clean-data/full-takeup-data.rds")
  

stop()



load(file.path("data", "analysis.RData"))
census.data

takeup_df = takeup.data %>%
  as_tibble()

# Loading census data
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

anne_census_df = census_data %>%
  filter(!is.na(cluster.id))

# 1352
anne_census_df %>%
  filter(true.monitored) %>%
  count(sms.treatment)
# 7155 + 2650 = 9805


anne_census_df %>%
  count(have_phone)

anne_census_df %>%
  filter(sms.treatment == "sms.control") %>%
  filter(monitored) %>%
  filter(have_phone == "No" | have_phone == "Don't know number") 

anne_census_df %>%
  filter(sms.treatment == "sms.control") %>%
  filter(true.monitored) %>%
  filter(have_phone == "No" | have_phone == "Don't know number") 


census.data %>%
  count(sms.treatment)


census_smsctrl_df = census.data %>%
  filter(is.na(sms.treatment) | sms.treatment == "sms.control") 

# 2,934
census_smsctrl_nophone_mon_df = census_smsctrl_df %>%
  filter(
    have_phone == "No" & monitored
  )

# 7,550
census_smsctrl_phone_mon_df = census_smsctrl_df %>%
  filter(have_phone == "Yes" & monitored)


census_smsctrl_df %>%
  filter(have_phone == "No") %>%
  count(monitored, true.monitored)
#
census.data %>%
  filter(is.na(sms.treatment) | sms.treatment == "sms.control")  %>%
  count(monitored, true.monitored)


mon_truemon_diff_keys = census_data %>%
  filter(is.na(sms.treatment) | sms.treatment == "sms.control")  %>%
  filter(monitored == TRUE & true.monitored == FALSE)  %>%
  pull(KEY.individ)




monitored_nosms_data %>%
  filter(KEY.individ %in% mon_truemon_diff_keys) 

census_data %>%
  summarize(n_distinct(cluster.id))




census_data %>%
  filter(is.na(sms.treatment) | sms.treatment == "sms.control")  %>%
  count(monitored, true.monitored)


census.data %>%
  filter(is.na(sms.treatment) | sms.treatment == "sms.control")  %>%
  filter(have_phone == "No") %>%
  filter(monitored == TRUE) %>%
  count(monitored, true.monitored)

census_data %>%
  filter(is.na(sms.treatment) | sms.treatment == "sms.control")  %>%
  filter(have_phone == "No") %>%
  filter(monitored == TRUE) %>%
  count(monitored, true.monitored)

stop()


  dewormed.day.data <- takeup.data %>% 
    filter(!is.na(KEY.individ)) %>% 
    rename(dewormed.day = deworming.day) %>% 
    group_by(KEY.individ) %>% 
    summarize(across(c(dewormed.day, dewormed.date), min)) %>% # If dewormed multiple times, take the first day only
    ungroup()
  
  rep_analysis.data <- census.data %>% 
    filter(!is.na(wave)) %>%  # Remove clusters no longer in study
    # left_join(.consent.dewormed.reports, by = "KEY.individ") %>% 
    mutate(
      dewormed = KEY.individ %in% discard(takeup.data$KEY.individ, is.na), # TRUE if individual found in take-up data
      dewormed = if_else(true.monitored, dewormed, NA), # We don't know the deworming status of those not monitored
      # hh.baseline.sample = KEY %in% .baseline.data$KEY, # We don't have a clear link between baseline surveys and the census. See the field notebook
      # monitor.consent = !is.na(monitor.consent) & monitor.consent,
      sms.treated = fct_match(sms.treatment, c("social.info", "reminder.only")),
      sms.treatment.2 = fct_explicit_na(sms.treatment, "sms.control") %>% fct_relevel("sms.control"),
      have.phone.bool = fct_match(have_phone, "Yes"),
      sms.ctrl.subpop = if_else(is.na(sms.ctrl.subpop) & !sms.treated, 
                                if_else(have.phone.bool, "phone.owner", "non.phone.owner"),
                                sms.ctrl.subpop)
    ) %>% # NA if not in the monitored group
    left_join(dewormed.day.data, by = "KEY.individ") %>% 
    filter(!sms.treated | (!is.na(sms.consent) & sms.consent)) # Drop those assigned to SMS treatment but were not consented (hence not treated)

  rep_analysis.data %>%
    filter(sms.treatment.2 == "sms.control") %>%
    count(true.monitored)  

  analysis.data %>% 
    filter(is.na(dewormed) | !dewormed) %>% # For anyone in study with with unknown or negative deworming status
    group_by(cluster.id) %>% 
    do(name.match.monitored(., filter(takeup.data, cluster.id %in% .$cluster.id), max.cost = max.name.match.cost)) %>% 
    ungroup %>% 
    right_join(analysis.data, by = c("cluster.id", "KEY.individ")) %>% 
    left_join(transmute(takeup.data, KEY.survey.individ, dewormed.day.matched = deworming.day, dewormed.date.matched = dewormed.date), 
              by = c("which.min.name.match.dist" = "KEY.survey.individ")) %>% 
    # mutate(monitored = !is.na(monitored) & !is.na(wave) & monitored, # Remove those dropped from the study 
    mutate(monitored = !is.na(true.monitored) & !is.na(wave) & true.monitored, # Remove those dropped from the study 
           dewormed.any = (!is.na(dewormed) & dewormed) | dewormed.matched,
           dewormed.day.any = if_else(!is.na(dewormed.day), as.integer(dewormed.day), dewormed.day.matched), 
           dewormed.date.any = if_else(!is.na(dewormed.date), dewormed.date, dewormed.date.matched),
           gender = factor(gender, levels = 1:2, labels = c("male", "female"))) %>% 
    left_join(
        transmute(.endline.data, 
                  KEY.individ, age, school, floor, ethnicity, ethnicity2, any.sms.reported, gift_choice, hh_cal, cal_value, hh_bracelet, number_bracelet, 
                  endline_deworm_rate = dworm_rate), by = "KEY.individ") %>% 
    mutate_at(vars(age, age.census), list(squared = ~ (.)^2, group = ~ cut(., breaks = c(seq(18, 58, 10), 120), right = FALSE))) %>% 
    left_join(select(.cluster.strat.data, wave, county, cluster.id, dist.pot.group), by = c("wave", "county", "cluster.id")) %>% 
    `attr<-`("class", c("takeup_df", class(.))) %>%
    unite(county_dist_stratum, county, dist.pot.group, remove = FALSE) %>% 
    unite(county_dist_mon_stratum, county, dist.pot.group, true.monitored, remove = FALSE) %>% 
    mutate(across(c(county, county_dist_stratum, county_dist_mon_stratum), factor)) %>% 
    mutate(mon_status = factor(true.monitored, levels = c(TRUE, FALSE), labels = c("monitored", "unmonitored")),
           name_matched = !true.monitored,
           phone_owner = !fct_match(sms.treatment.2, "sms.control") | fct_match(sms.ctrl.subpop, "phone.owner"))


  1352 + 5803


stop()

