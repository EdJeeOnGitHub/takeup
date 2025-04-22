library(tidyverse)
library(here)

baseline.data = read_rds("temp-data/reclean_baseline_data.rds") # Not sampling data!




baseline_data = read_rds(file.path("data", "takeup_baseline_data.rds"))

all_endline_data = read_rds(file.path("data", "all_endline.rds"))

baseline_data



all_endline_data %>%  
  select(religion, ethnicity)

all_endline_data %>%
  count(ethnicity) %>%
  arrange(-n)

endline_data %>%
  select(contains("clan"))

baseline_data %>%
  select(religion, ethnicity) %>%
  unique()

all_endline_data

# Baseline + Endline
# merge endline:
  # - phone
  # - primary
  # - floor

  # top two from each
  # - ethnicity: proportion majority ethnicity
  # - religion: proportion majority religion



# Clean up know everyone can be infected
# knowledge variables + add more sample

all_endline_data %>%
  select(contains('floor'))

all_endline_data %>%
  select(contains('phone'))

all_endline_data %>%
  select(contains('school')) %>%
  count(school) %>%
  arrange(-n)




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

census_data %>%
  colnames()

    census_data %>%
      summarise(mean(have_phone == "Yes", na.rm = TRUE))

stop()
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

stop()

raw_baseline_data = raw_baseline_data %>%
  mutate(
    sub_date_filter = SubmissionDate >= "2016-09-05",
    present_filter = present == 1 | !is.na(age),
    consent_filter = !is.na(consent) & consent == 1
  ) 
  %>%
  summarise(
    across(contains('filter'), mean, na.rm = TRUE)
  )


raw_baseline_data %>%
  count(sub_date_filter)
# GOAL 2160
# HAVE 2056

raw_baseline_data %>%
  filter(sub_date_filter == TRUE) %>%
  count(present_filter) 

raw_baseline_data %>%
  filter(sub_date_filter == TRUE) %>%
  filter(present_filter == TRUE) %>%
  count(consent_filter)

  filter(
    sub_date_filter
  )

  filter(SubmissionDate >= "2016-09-05",
         present == 1 | !is.na(age), 
         !is.na(consent) & consent == 1) 

reclean_baseline_data = tu.data.reader(raw.data.path("Baseline Survey.csv")) %>% 
  filter(SubmissionDate >= "2016-09-05", 
         present == 1 | !is.na(age), 
         !is.na(consent) & consent == 1) %>% 
  select(-county)  %>%
  left_join(filter(., !invalid.coord) %>% identify.closest.cluster, "KEY") %>%
  left_join(cluster.wave.county.data, "cluster.id") 


