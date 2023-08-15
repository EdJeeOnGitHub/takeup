library(magrittr)
library(tidyverse)
library(lubridate)
library(haven)
library(broom)
library(scales)
library(rgeos)
library(here)


source('analysis_util.R')
source("rct-design-fieldwork/takeup_rct_assign_clusters.R")


wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

datetime.format <- "%b %d, %Y %I:%M:%S %p"
takeup.datetime.type <- col_datetime(datetime.format)
takeup.date.type <- col_date(datetime.format)
raw.data.path <- . %>% here("data", "raw-data", .)

rct.counties <- c("Busia", "Siaya", "Kakamega")
busia.subcounties <- c("butula", "nambale", "teso south", "teso north") 
siaya.subcounties <- c("gem", "ugenya", "ugunja")

ke.lvl2.adm.data <- read_rds(here("data", "adm", "KEN_adm2.rds"))
#ke.lvl3.adm.data <- read_rds("~/Data/TakeUp/KEN_adm3.rds")

counties.adm.data <- ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ] #, "Vihiga"), ]
subcounties.adm.data <- counties.adm.data[!counties.adm.data$NAME_1 %in% c("Busia", "Siaya") | counties.adm.data$NAME_2 %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ]


clusters.to.drop <- c(277, # Too close other cluster (I think 503)
                      491, 492, # Problematic urban clusters
                      1, # Village dispute about PoT
                      678, # Hostile community member
                      737) # Data fabrication and medication theft

cluster.wave.county.data <- read_rds(here("data", "takeup_cluster_wave_county_5.0.rds"))

cluster.info.to.drop <- c("uuid:3a491628-a0e8-4be0-9afc-c2cca08fc450", # Bad entries
                          "uuid:fcd3137e-7ddb-4f2e-9140-680e7638e42c",
                          "uuid:00fe16c6-5de1-4ff8-a9c8-92101ef4bd01",
                          "uuid:e3a47c54-f7f2-4af0-bd5b-0de12d997cb0")

cluster.info <- read_csv(raw.data.path("Cluster Survey V3 July 04.csv"),
                         col_types = list(SubmissionDate = takeup.date.type)) %>% 
  bind_rows(v3 = ., v1 = read_dta(here("data", "Cluster Survey V1.dta")) %>%
                   transmute(cluster.id = clusterid,
                             SubmissionDate = submissiondate,
                             KEY = key,
                             `gps2-Longitude` = gps2longitude,
                             `gps2-Latitude` = gps2latitude,
                             SubmissionDate = parse_date(SubmissionDate, "%d/%m/%Y %T"),
                             deviceid = as.character(deviceid),
                             manual_long2, manual_lat2, location_type, alt_name, comments, enumerator), 
            .id = "data.source") %>%
  filter(!is.na(cluster.id), deviceid != "(web)", !KEY %in% cluster.info.to.drop) %>% 
  arrange(cluster.id, SubmissionDate) %>% 
  group_by(cluster.id) %>%
  slice(n()) %>%
  ungroup %>%
  rename(alt.pot.lon = `gps2-Longitude`,
         alt.pot.lat = `gps2-Latitude`) %>%
  mutate(alt.pot.lon = ifelse(is.na(alt.pot.lon), read_lon_ant, alt.pot.lon),
         alt.pot.lat = ifelse(is.na(alt.pot.lat), read_lat_ant, alt.pot.lat),
         alt.pot.lon = ifelse(is.na(alt.pot.lon), manual_long2, alt.pot.lon),
         alt.pot.lat = ifelse(is.na(alt.pot.lat), manual_lat2, alt.pot.lat),
         location_type = factor(location_type, levels = 1:5, labels = c("Clinic", "Church", "Market", "Home", "Other")),
         cluster.id = as.integer(cluster.id)) %>% 
  rename(old.county.code = county) %>%
  mutate(old.county.code = factor(old.county.code, levels = 1:3, labels = c("Busia", "Kakamega", "Siaya"))) %>% 
  left_join(cluster.wave.county.data, "cluster.id") 

pot.info <- cluster.info %>%
  select(enumerator, enumerator_other, wave, county, cluster.id, alt.pot.lon, alt.pot.lat, alt_name, location_type, comments, SubmissionDate, data.source)

pot.verify.data <- read_csv(raw.data.path("POT verification.csv"), 
                            col_types = list(SubmissionDate = takeup.datetime.type)) %>%
  rename(cluster.id = cluster_id,
         lon.verify = `gps-Longitude`,
         lat.verify = `gps-Latitude`) %>% 
  select(-county) %>% 
  left_join(cluster.wave.county.data, "cluster.id") %>% 
  filter(wave == 1 | SubmissionDate >= "2016-10-16", !is.na(lon.verify), !is.na(lat.verify)) %>% 
  group_by(cluster.id) %>% 
  mutate(num.entries = n()) %>% 
  filter(min_rank(SubmissionDate) == n()) %>% 
  ungroup 

pot.info %<>% 
  left_join(#filter(pot.verify.data, !is.na(lon.verify), !is.na(lat.verify)),
            pot.verify.data,
            c("wave", "county", "cluster.id"), 
            suffix = c(".original", ".verify")) %>% 
  set_names(str_replace(names(.), "\\.original$", "")) 

pot.info %<>% 
  filter(!is.na(lon.verify), !is.na(lat.verify), !is.na(alt.pot.lon), !is.na(alt.pot.lat)) %>% 
  mutate(verify.dist = rgeos::gDistance(convert.to.sp(., ~ alt.pot.lon + alt.pot.lat, wgs.84) %>% spTransform(kenya.proj4),
                                 convert.to.sp(., ~ lon.verify + lat.verify, wgs.84) %>% spTransform(kenya.proj4), byid = TRUE) %>% diag) %>% 
  select(KEY, verify.dist) %>% 
  right_join(pot.info, "KEY")







county.bbox <- counties.adm.data@bbox
subcounty.bbox <- subcounties.adm.data@bbox
census.data = read_rds("data/takeup_census.rds")

census_dict <- census.data %>% 
  mutate_at(vars(name1st, name_mid, name2nd, hhh_name2nd, hhh_full_name), . %>% coalesce("") %>% str_trim() %>% str_to_upper()) %>% 
  transmute(
    cluster.id, village, KEY, KEY.individ, hhh_name1st, hhh_name2nd, name1st, age.census,
    name_census = sprintf("%s %s %s", name1st, name_mid, name2nd) %>% str_replace_all("\\s+", " "),
    hhh_name = hhh_full_name %>% str_replace_all("\\s+", " "))


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

rct.village.codes <- read_csv(here("data", "village_codes_2.csv"), skip = 1, col_names = c("village.cluster.id", "village_name", "village")) %>% 
  mutate(village.cluster.id = as.integer(village.cluster.id))

# We're already reading this above (cluster-survey-data), but not in exactly the same format...
cluster.survey.data <- read_rds(here("data", "takeup_cluster_survey.rds"))

rct.villages <- read_rds(here("data", "rct_target_villages_2.0.rds")) %>% 
  mutate(new.village = FALSE) %>% 
  bind_rows(read_rds("data/rct_target_villages_2.0-4.rds") %>% 
              mutate(new.village = TRUE))
  
all.villages <- rct.villages %>% 
  bind_rows(anti_join(cluster.survey.data, ., c("cluster.id", "target.village.id"))) %>% 
  mutate(targeted.village = !is.na(new.village),
         cluster.id = as.integer(cluster.id),
         target.village_name = str_trim(target.village_name) %>% str_replace_all("\\s+", " ")) %>% 
  left_join(rct.village.codes, c("cluster.id" = "village.cluster.id", "target.village_name" = "village_name")) %>% 
  left_join(select(pot.info, cluster.id, matches("alt.pot.(lon|lat)"), location_type), "cluster.id") %>% 
  group_by(is.na(village)) %>% 
  do({ # Adding some village IDs for village that weren't in the code file shared by Arthur
    if(is.na(first(.$village))) {
      mutate(., village = 1000 + seq_len(nrow(.)))
    } else {
      (.)
    }
  }) %>% 
  ungroup %>% 
  mutate(dist.group = convert.to.sp(., ~ target.lon + target.lat, wgs.84) %>% 
           spTransform(kenya.proj4) %>% 
           gDistance(byid = TRUE) %>% 
           as.dist %>% 
           hclust %>% 
           cutree(h = 1000)) %>% 
  group_by(cluster.id) %>% 
  mutate(num.dist.groups = unique(dist.group) %>% length) %>% 
  ungroup %>% 
  left_join(cluster.wave.county.data, "cluster.id") %>% 
  mutate(village.name.group = target.village_name %>% 
           str_trim %>% 
           str_replace("\\s*[A-Z]\\d?$", "") %>% 
           str_replace_all("\\s+", " ") %>% 
           str_replace_all("'", "") %>% 
           str_to_upper %>% 
           str_replace("\\s+[A-Z]\\d?$", "") %>% 
           str_replace(regex("\\s+(village|upper|lower|group|east|west|north|south|rural|urban|township|central)$", ignore_case = TRUE), "") %>% 
           str_replace(regex("\\s+(MWILUECHINA)$", ignore_case = TRUE), "") %>% 
           str_replace(regex("\\s*(\\(.+\\)|estate|township)$", ignore_case = TRUE), "") %>% 
           str_replace("\\s+[A-Z]\\d?$", ""),
         vill.name.dist.group = adist(village.name.group) %>% 
           as.dist %>% 
           hclust %>% 
           cutree(h = 1))

rm(rct.villages)

rct.village.codes %<>%
  left_join(select(all.villages, village, targeted.village), "village")



known.pot.locations <- pot.info %>%
  filter(!is.na(alt.pot.lon) | !is.na(lon.verify), !is.na(alt.pot.lat) | !is.na(lat.verify)) %>%
  mutate(alt.pot.lon = ifelse(is.na(alt.pot.lon), lon.verify, alt.pot.lon),
         alt.pot.lat = ifelse(is.na(alt.pot.lat), lat.verify, alt.pot.lat)) %>% 
  distinct(cluster.id, .keep_all = TRUE)

known.village.locations <- all.villages %>%  #rct.villages %>%
  filter(targeted.village, !is.na(target.lon), !is.na(target.lat)) %>%
  distinct(cluster.id, target.village.id, .keep_all = TRUE)

identify.closest.cluster <- function(.data, data.coords.formula = ~ lon + lat,  key.variable = "KEY") {
  # min.dist.df <- .data %>% 
  .data %>%
    filter(rowSums(is.na(model.frame(data.coords.formula, data = ., na.action = NULL))) == 0) %>% 
    # filter(!is.na(lon), !is.na(lat)) %>%
    (function (.located.data) {
      convert.to.sp(.located.data, data.coords.formula, wgs.84) %>% 
        spTransform(kenya.proj4) %>% 
        gDistance(known.village.locations %>%
                    convert.to.sp(~ target.lon + target.lat, wgs.84) %>%
                    spTransform(kenya.proj4),
                  byid = TRUE) %>%
        plyr::adply(2, function(dist.col) {
          in.range.dist <- dist.col[dist.col <= 2000]
      
          if (length(in.range.dist) == 0) {
            tibble(min.dist = NA, closest.cluster = NA)
          } else {
            tibble(min.dist = min(in.range.dist),
                   closest.village = known.village.locations$target.village.id[which(dist.col <= min.dist)], 
                   closest.cluster = known.village.locations$cluster.id[which(dist.col <= min.dist)]) %>% 
              distinct(closest.cluster, .keep_all = TRUE)
          }
          
        }, .parallel = TRUE) %>% 
        bind_cols(.located.data, .) %>% 
        select_(.dots = c(key.variable, "min.dist", "closest.cluster")) 
    })
}


raw_baseline_data <- tu.data.reader(raw.data.path("Baseline Survey.csv")) 


reclean_baseline_data = tu.data.reader(raw.data.path("Baseline Survey.csv")) %>% 
  filter(SubmissionDate >= "2016-09-05", 
         present == 1 | !is.na(age), 
         !is.na(consent) & consent == 1) %>% 
  select(-county)  %>%
  left_join(filter(., !invalid.coord) %>% identify.closest.cluster, "KEY") %>%
  left_join(cluster.wave.county.data, "cluster.id") 




baseline.data <- tu.data.reader(raw.data.path("Baseline Survey.csv")) %>% 
  filter(SubmissionDate >= "2016-09-05", 
         !is.na(present) & present == 1, 
         !is.na(consent) & consent == 1) %>% 
  select(-county) %>% 
  left_join(filter(., !invalid.coord) %>% identify.closest.cluster, "KEY") %>%
  left_join(cluster.wave.county.data, "cluster.id") %>% 
  filter(!is.na(wave)) 

reclean_baseline_data %<>% 
  mutate(name_census = str_replace_all(name_census, "\\s+", " ") %>% str_trim(),
         hhh_name = str_replace_all(hhh_name, "\\s+", " ") %>% str_to_upper() %>% str_trim(),
         hhh_first_name = str_extract(hhh_name, "^\\w+"),
         hhh_last_name = str_extract(hhh_name, "\\S+$")) 


cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds")) # Study clusters metadata
reclean_baseline_data %<>% 
  prepare.baseline.data(cluster.strat.data) 

write_rds(
  reclean_baseline_data, 
  "temp-data/reclean_baseline_data.rds") # Not sampling data!



reclean_baseline_data
