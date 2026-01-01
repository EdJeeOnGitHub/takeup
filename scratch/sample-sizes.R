library(tidyverse)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

# Data --------------------------------------------------------------------

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

# stick to monitored sms.treatment group
# remove sms.treatment.2

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


monitored_sms_data = analysis.data %>% 
  filter(mon_status == "monitored") %>% 
  filter(have_phone == "Yes") %>%
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

baseline.data %>%
  write_csv("temp-data/debug-baseline-data.csv")


endline.data %>%
  write_csv("temp-data/debug-endline-data.csv")


endline.data %>%
  group_by(sms.treatment) %>%
  summarise(
    n_hh = n_distinct(KEY)
  )






monitored_sms_data


monitored_nosms_data

analysis.data

analysis.data %>%
  count(sms.treatment.2, monitored)

analysis.data %>%
  count(sms.treatment.2)


monitored_sms_data %>%
  count(sms.treatment)


monitored_sms_data %>%
  count(sms.treatment.2) 

monitored_sms_data %>%
  count(sms.treatment, sms.treatment.2)


#### SMS Stuff
monitored_nosms_data

monitored_sms_data %>%
    count(sms.treatment)

monitored_data = analysis.data %>%
    filter(mon_status == "monitored") 

monitored_data %>%
  count(wave)


monitored_data %>%
    select(contains("consent")) %>%
    skimr::skim()


nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


analysis_data <- monitored_nosms_data %>% 
  ungroup() %>%
  mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group,
    standard_hh.dist.to.pot = standardize(dist.to.pot),
    indiv_id = row_number()
    )

### WTP Number of observations
wtp_stan_data <- analysis.data %>% 
  mutate(stratum = county) %>% 
  prepare_bayes_wtp_data(
    wtp.data,
    
    preference_value_diff = seq(-100, 100, 10), 
    num_preference_value_diff = length(preference_value_diff), 
    
    wtp_utility_df = 3,
    tau_mu_wtp_diff = 100,
    mu_wtp_df_student_t = 7,
    tau_sigma_wtp_diff = 50,
    sigma_wtp_df_student_t = 2.5
  )

wtp_out = analysis.data %>%
  mutate(stratum = county) %>%
  prepare_bayes_wtp_data(
    wtp.data,
    
    preference_value_diff = seq(-100, 100, 10), 
    num_preference_value_diff = length(preference_value_diff), 
    
    wtp_utility_df = 3,
    tau_mu_wtp_diff = 100,
    mu_wtp_df_student_t = 7,
    tau_sigma_wtp_diff = 50,
    sigma_wtp_df_student_t = 2.5
  )

  wtp.data  %>%
    summarise(
      n = n_distinct(KEY)
    ) 

analysis.data %>%
    filter(
           assigned.treatment == "control") %>%
    count(
      wtp_samp = !is.na(gift_choice), 
      sms_control = sms.treatment.2 == "sms.control",
      monitored = true.monitored
      ) 

      wtp.data

analysis.data %>%
    filter(
           assigned.treatment == "control",
           sms.treatment.2 == "sms.control")  %>%
    count(gift_choice) %>%
    filter(!is.na(gift_choice)) %>%
    summarise(
      n = sum(n)
    )



# Gift preference for bracelet vs calendar sample
analysis_data %>%
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
    )  


analysis_data %>%
  # 9805
  filter(!is.na(gift_choice)) %>%
  # 2,312
  filter(monitored) %>%
  # 2,312
  filter(monitor.consent)  %>%
  # 2,312
  filter(!hh.baseline.sample.pool)  %>%
  # 1940
  filter(!is.na(sms.treatment))  %>%
  # 1940
  filter(dewormed == FALSE)
  # 1,174


# Endline ------------------------------------------------------------------


# total monitored
analysis.data %>%
  filter(mon_status == "monitored")

# total monitored w/ SMS Control
analysis.data %>%
  filter(mon_status == "monitored")  %>%
  filter(sms.treatment.2 == "sms.control")
# total monitored in SMS
analysis.data %>%
  filter(mon_status == "monitored")  %>%
  filter(sms.treatment.2 != "sms.control")



total = 9935
baseline = 2250
endline = 3750
sms = 3935

all.endline.data %>%
  count(consent)

all.endline.data %>%
  select(KEY.individ)

baseline.data %>%
  select(KEY) %>%
  slice(1:3)

monitored_sms_data %>%
  select(KEY) %>%
  slice(1:3)

baseline.data %>%
  colnames()

intersect(
  baseline.data %>%
    pull(KEY),
  census.data %>%
    pull(KEY)
)

baseline.data


census.data %>%
  filter(true.monitored == 1) %>%
  group_by(cluster.id, sms.treatment, assigned.treatment) %>%
  summarise(
    n_hh = n_distinct(KEY),
    n_indiv = n_distinct(KEY.individ)
  )

monitored_nosms_data %>%
  count(hh.baseline.sample.pool)

monitored_sms_data %>%
  count(hh.baseline.sample.pool)

analysis_data %>%
  count(hh.baseline.sample.pool)

analysis_data %>%
    filter(!is.na(gift_choice), monitored, monitor.consent, !hh.baseline.sample.pool, !is.na(sms.treatment)) 

monitored_sms_data %>%
  select(KEY.individ)


baseline + endline



baseline + endline + sms


baseline.data

endline.know.table.data %>%
    filter(fct_match(know.table.type, "table.A"))  %>%
    group_by(KEY.individ) %>%
    mutate(num.recognized = sum(num.recognized))  %>%
    filter(num.recognized > 0) %>%
    ungroup() %>%
    summarise(
      n_indiv = n_distinct(KEY.individ)
    )


know_a_table = endline.know.table.data %>%
  filter(fct_match(know.table.type, "table.A")) 

sms_treat = monitored_sms_data %>%
  filter(sms.treatment %in% c("social.info", "reminder.only")) 

sms_treat


sms_treat

# out of total endline sample, we randomly sampled 1627 people to implement beliefs.
# Amongst them, XYZ didn't recognize anyone, and remaining ABC in the SMS sample, 
# who are not used in the main analysis.o


# out of endline sample, we randomly sampled 1627 people to implement beliefs.
# In the non-SMS sample this leads to ABC total beliefs observations, of which 
# XYZ didn't recognize anyone.


know_a_clusters = know_a_table %>%
  pull(cluster.id)

ana_clusters = analysis_data %>%
  pull(cluster.id)

setdiff(ana_clusters, know_a_clusters)


know_a_ids = know_a_table %>%
  pull(KEY.individ)

ana_ids = analysis_data %>%
  pull(KEY.individ)

setdiff( know_a_ids, ana_ids)

not_in_monitored = setdiff(know_a_ids, monitored_data$KEY.individ)

not_in_monitored %>%
  saveRDS("temp-data/not_in_monitored.rds")


all.endline.data %>%
    filter(across(c(present, interview, consent), ~ !is.na(.x) & .x == 1)) 

# interviewied 3817, of which 3716 consented
all.endline.data %>%
    filter(present == TRUE) %>%
    count(consent) 


# of which 3700 finished interview
all.endline.data %>%
    filter(present == TRUE & consent == TRUE) %>%
    count(interview)

# of which 3678 weren't duplicate people?!
all.endline.data %>% 
    filter(present == TRUE & consent == TRUE & interview == TRUE) %>%
    group_by(KEY.individ) %>%
    filter(row_number() == 1)  %>%
    nrow()


all.endline.data %>%
    filter(across(c(present, interview, consent), ~ !is.na(.x) & .x == 1)) %>% 
    arrange(KEY.individ, SubmissionDate) %>% 
    group_by(KEY.individ) %>% 
    filter(row_number() == 1)  # If more than one entry for an individual, take first one (there are 22 such individuals)

stop()

all.endline.data %>%
  count(consent)

  colnames(all.endline.data)


all.endline.data %>%
  select(present, return, consent, reconsent)

all.endline.data %>%
  count(reconsent)


endline.data

know_a_table = know_a_table %>%
  group_by(KEY.individ) %>%
  mutate(num.recognized = sum(num.recognized)) %>%
  mutate(any_recognized = num.recognized > 0) %>%
  ungroup() %>%
  mutate(
    id_in_analysis = if_else(KEY.individ %in% analysis_data$KEY.individ, 1, 0),
    id_in_sms = if_else(KEY.individ %in% sms_treat$KEY.individ, 1, 0),
    id_in_either = if_else(KEY.individ %in% analysis_data$KEY.individ | KEY.individ %in% sms_treat$KEY.individ, 1, 0),
    id_in_mon = if_else(KEY.individ %in% monitored_data$KEY.individ, 1, 0)
  )

know_a_table %>%
  summarise(
    pr_in_ana = mean(id_in_analysis, na.rm = TRUE),
    pr_in_sms = mean(id_in_sms, na.rm = TRUE),
    pr_in_either = mean(id_in_either, na.rm = TRUE)
  )

know_a_table %>%
  group_by(id_in_analysis, id_in_sms, id_in_mon) %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )

in_know_not_in_mon = know_a_table %>%
  filter(id_in_mon == 0) %>%
  pull(KEY.individ)

endline.data %>%
  filter(
    KEY.individ %in% in_know_not_in_mon
  ) 
endline.data %>%
  filter(
    KEY.individ %in% in_know_not_in_mon
  ) %>%
  group_by(assigned.treatment) %>%
  summarise(n = n()) %>%
  count()

all.endline.data %>%
  filter(
    KEY.individ %in% in_know_not_in_mon
  )



know_a_table %>%
  filter(id_in_mon == 0)  %>%
  count(dewormed)


know_a_table %>%
  filter(id_in_mon == 1)  %>%
  count(dewormed)

know_a_table %>%
  filter(id_in_analysis == 1) %>%
  group_by(KEY.individ) %>%
  summarise(
    n_recognized = sum(num.recognized)
  ) %>%
  count(n_recognized) 

know_a_table %>%
  group_by(id_in_analysis, id_in_sms, any_recognized) %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )



#### Baseline Data 

library(tidyverse)
library(here)
library(sp)

baseline.data = read_rds("temp-data/reclean_baseline_data.rds") # Not sampling data!
baseline_data = read_rds(file.path("data", "takeup_baseline_data.rds"))
all_endline_data = read_rds(file.path("data", "all_endline.rds"))




# Clean up know everyone can be infected
# knowledge variables + add more sample

# 9,935 adults of which 2,250 adults are surveyed at baseline, 
# 3,750 adults surveyed at endline and 
# 3,935 adults selected for text messaging intervention.


print(str_glue("Endline data: {nrow(all_endline_data)}"))
print(str_glue("Baseline data: {nrow(baseline_data)}"))
nrow(baseline.data)

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


n_hh_df = census_data %>%
    group_by(cluster.id) %>%
    summarise(
        n = sum(num.individuals)
    )

rct.counties <- c("Busia", "Siaya", "Kakamega")
busia.subcounties <- c("butula", "nambale", "teso south", "teso north") 
siaya.subcounties <- c("gem", "ugenya", "ugunja")

ke.lvl2.adm.data <- read_rds(here("data", "adm", "KEN_adm2.rds"))
#ke.lvl3.adm.data <- read_rds("~/Data/TakeUp/KEN_adm3.rds")

counties.adm.data <- ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ] #, "Vihiga"), ]
subcounties.adm.data <- counties.adm.data[!counties.adm.data$NAME_1 %in% c("Busia", "Siaya") | counties.adm.data$NAME_2 %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ]

datetime.format <- "%b %d, %Y %I:%M:%S %p"
takeup.datetime.type <- col_datetime(datetime.format)
takeup.date.type <- col_date(datetime.format)
raw.data.path <- . %>% here("data", "raw-data", .)
#### Baseline data ####
validate.coords <- . %>% 
  mutate(invalid.coord = 
           (!is.na(lon) & (lon > county.bbox["x", "max"] | lon < county.bbox["x", "min"])) |
           (!is.na(lat) & (lat > county.bbox["y", "max"] | lat < county.bbox["y", "min"])))
tu.data.reader <- function(file.name, submit.datetime.type = NULL, .other.types = NULL) { # =  "%b %d, %Y %I:%M:%S %p") {
  col.types <- list(SubmissionDate = if (is.null(submit.datetime.type)) takeup.datetime.type else submit.datetime.type,
                                       manual_long = col_number(),
                                       manual_lat = col_number()) %>% 
    c(.other.types)
  
  read_csv(file.name, col_types = col.types) %>% 
    mutate(isValidated = isValidated == "true") %>% 
    rename(lat = `gps-Latitude`,
           lon = `gps-Longitude`,
           cluster.id = cluster_id) %>% 
    filter(deviceid != "(web)") %>% 
    mutate(lon = ifelse(is.na(lon), manual_long, lon),
           lat = ifelse(is.na(lat), manual_lat, lat)) %>% 
    validate.coords()
}
county.bbox <- counties.adm.data@bbox
subcounty.bbox <- subcounties.adm.data@bbox
census.data = read_rds("data/takeup_census.rds")

raw_baseline_data <- tu.data.reader(raw.data.path("Baseline Survey.csv")) 


raw_baseline_data = raw_baseline_data %>%
  mutate(
    sub_date_filter = SubmissionDate >= "2016-09-05",
    present_filter = present == 1 | !is.na(age),
    consent_filter = !is.na(consent) & consent == 1
  ) 

raw_baseline_data = raw_baseline_data %>%
  mutate(
    date = sub_date_filter,
    date_and_present = sub_date_filter & present_filter,
    date_and_present_and_consent = date_and_present & consent_filter
  )


# Baseline breakdown
# Actually talk to 3001 people - some of these seem to be too early - pilot?
raw_baseline_data %>%
  group_by(date) %>%
  summarise(
    n = n()
  )
# of which, only 2151 (goal was 2160) are present
raw_baseline_data %>%
  group_by(date_and_present) %>%
  summarise(
    n = n()
  )
# of which, only 2069 consent
raw_baseline_data %>%
  group_by(date_and_present_and_consent) %>%
  summarise(
    n = n()
  )

census.data %>%
  count(monitored, true.monitored)

census.data %>%
  count(true.monitored)


census.data %>%
  filter(true.monitored == 1) %>%
  count(recruit)

