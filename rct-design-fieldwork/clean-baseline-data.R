library(magrittr)
library(tidyverse)
library(sf)
library(here)

source("scripts/shared/clean-analysis-setup.R")
source("rct-design-fieldwork/takeup_rct_assign_clusters.R")
dir.create("data/clean-data", showWarnings = FALSE)

wgs.84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 = "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

datetime.format = "%b %d, %Y %I:%M:%S %p"
takeup.datetime.type = col_datetime(datetime.format)
raw.data.path = . %>% here("data", "raw-data", .)

rct.counties = c("Busia", "Siaya", "Kakamega")

ke.lvl2.adm.data = read_rds(here("data", "adm", "KEN_adm2.rds"))
counties.adm.data = ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ]

cluster.wave.county.data <- read_rds(here("data", "takeup_cluster_wave_county_5.0.rds"))

county.bbox <- counties.adm.data@bbox

validate.coords <- . %>%
  mutate(invalid.coord =
           (!is.na(lon) & (lon > county.bbox["x", "max"] | lon < county.bbox["x", "min"])) |
           (!is.na(lat) & (lat > county.bbox["y", "max"] | lat < county.bbox["y", "min"])))
cluster.survey.data <- read_rds(here("data", "takeup_cluster_survey.rds"))

rct.villages <- read_rds(here("data", "rct_target_villages_2.0.rds")) %>%
  mutate(new.village = FALSE) %>%
  bind_rows(read_rds("data/rct_target_villages_2.0-4.rds") %>%
              mutate(new.village = TRUE))

all.villages <- rct.villages %>%
  bind_rows(anti_join(cluster.survey.data, ., c("cluster.id", "target.village.id"))) %>%
  mutate(targeted.village = !is.na(new.village),
         cluster.id = as.integer(cluster.id))

rm(rct.villages)

known.village.locations <- all.villages %>%
  filter(targeted.village, !is.na(target.lon), !is.na(target.lat)) %>%
  distinct(cluster.id, target.village.id, .keep_all = TRUE)


reclean_baseline_data = tu_data_reader(raw.data.path("Baseline Survey.csv")) %>%
  filter(SubmissionDate >= "2016-09-05",
         present == 1 | !is.na(age),
         !is.na(consent) & consent == 1) %>%
  select(-county)  %>%
  left_join(filter(., !invalid.coord) %>% identify_closest_cluster, "KEY") %>%
  left_join(cluster.wave.county.data, "cluster.id")

reclean_baseline_data = reclean_baseline_data %>%
  mutate(name_census = str_replace_all(name_census, "\\s+", " ") %>% str_trim(),
         hhh_name = str_replace_all(hhh_name, "\\s+", " ") %>% str_to_upper() %>% str_trim(),
         hhh_first_name = str_extract(hhh_name, "^\\w+"),
         hhh_last_name = str_extract(hhh_name, "\\S+$"))

cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))

reclean_baseline_data = reclean_baseline_data %>%
  prepare_baseline_data(cluster.strat.data)


n_reached = reclean_baseline_data %>% filter(present == 1) %>% nrow()

n_reached = reclean_baseline_data %>%
  summarize(n = n_distinct(KEY)) %>%
  pull(n)

cat("Reached ", n_reached, " individuals in baseline survey\n")

# subset to those completed (wave not NA)
reclean_baseline_data = reclean_baseline_data %>%
  filter(!is.na(wave))

n_complete = reclean_baseline_data %>%
  summarize(n = n_distinct(KEY)) %>%
  pull(n)

cat("Complete data on ", n_complete, " individuals in baseline survey\n")

# 1995 observations
write_csv(
  reclean_baseline_data,
  "data/clean-data/clean-baseline-data.csv")

write_rds(
  reclean_baseline_data,
  "data/clean-data/clean-baseline-data.rds")
