# Load libraries ----------------------------------------------------------
library(plyr)       # for ldply - MUST load before tidyverse
library(tidyverse)
library(magrittr)
library(sf)
library(haven)      # for read_dta
library(lubridate)
raw_data_path <- "data"

# Coordinate reference systems --------------------------------------------
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

# Clean data --------------------------------------------------------------

# Load required data files
rct.targetable.schools <- read_rds("data/rct_targetable_schools_2.0.rds")

rct.schools.data <- read_rds(file.path(raw_data_path, "takeup_rct_schools.rds"))
if (!inherits(rct.schools.data, "sf")) {
  rct.schools.data <- st_as_sf(rct.schools.data)
}
rct.schools.data %<>% left_join(rct.targetable.schools, by = c("cluster.id" = "targeted.cluster.id"))

# Read and convert to sf if needed
rct.cluster.selection <- read_rds(file.path(raw_data_path, "rct_cluster_selection_2.0.rds"))
if (!inherits(rct.cluster.selection, "sf")) {
  rct.cluster.selection <- st_as_sf(rct.cluster.selection)
}

# Extract coordinates and join with school data
rct.cluster.selection <- rct.cluster.selection %>%
  left_join(rct.schools.data %>%
              st_transform(wgs.84) %>%
              mutate(lon = st_coordinates(.)[, "X"],
                     lat = st_coordinates(.)[, "Y"]) %>%
              st_drop_geometry() %>%
              select(cluster.id, lon, lat)) %>%
  rename(pot.lon = lon, pot.lat = lat)

read_cluster_survey_data <- function(data.filename) {
    df = read_dta(data.filename) %>% {
    if (!any(grepl("^anchor_gps", names(.)))) mutate(., anchor_gpslatitude = NA, anchor_gpslongitude = NA) else .
  } %>% {
    if (!"correct_ward" %in% names(.)) mutate(., correct_ward = NA, correction = NA) else .
  } %>% #{
    rename(cluster.id = clusterid) %>% 
    select(cluster.id, found, surveydate, new_name, manual_lat, manual_long,
           correct_ward, correction, # Ward updates
           boda, boda_contact, alt_boda, alt_boda_contact,
           matches("^gps[2-3]?\\.?l(at|on)"), starts_with("gps_work"),
           alt_name, location_type, other_type,
           # starts_with("market_days"),
           starts_with("anchor_gps"),
           matches("^(village_name|(second_)?(elder|chv)\\d?(_(contact|phone2|2nd_phone|2nd_contact|name))?|manual_.+|(?#anchor_(big|local)_market_[1-7]|)pop_households|pop_individuals)_?\\d"),
           matches("^pop_(households|individuals)_?\\d")) %>% 
    rename(anchor.gps.work = gps_work2,
           anchor.manual.lat = manual_lat2,
           anchor.manual.lon = manual_long2) %>%
    set_names(sub("^anchor_(.+\\d)$", "target.\\1", names(.))) %>% 
    set_names(sub("(.+)(\\d)(\\d)$", "\\1\\2_\\3", names(.))) %>% 
    set_names(sub("(village_name|gps3.+[^_]|check|pop_.*[^_])_?(\\d)$", "\\1_\\2", names(.))) %>% 
    set_names(ifelse(grepl("^(?!market|target).+_?\\d{1,2}$", names(.), perl = TRUE), paste0("target.", names(.)), names(.))) %>% 
    gather(key, val, starts_with("target.")) %>% 
    tidyr::extract(key, c("var", "target.village.id"), "(.+?)_?(\\d+)$") %>% 
    spread(var, val) %>% 
    filter(!is.na(target.village_name), nchar(target.village_name) > 0) %>% 
    mutate_at(vars(matches("target.+(long|lat|pop_(households|individuals))")), ~as.numeric(.)) %>% 
    mutate(data.filename = data.filename,
        #    surveydate = lubridate::as_date(surveydate),
        surveydate = as.character(surveydate),
           corrected.ward = ifelse(correct_ward == 1, correction, NA),
           corrected.pot.school.name = ifelse(nchar(new_name) > 0, new_name, NA)) %>% 
    select(-c(correct_ward, correction, new_name)) %>% 
    rename(alt.pot.name = alt_name,
           alt.pot.location.type = location_type,
           alt.pot.location.type.other = other_type)  %>%
    mutate(
        manual_lat = as.character(manual_lat),
        manual_long = as.character(manual_long),
        boda = as.character(boda),
        alt_boda = as.character(alt_boda)
    )
    
    return(df)
}


df = map_dfr(c("Cluster Survey V1.dta", "Cluster Survey 27.07.2016.dta"), ~read_cluster_survey_data(file.path(raw_data_path, .))) 


cluster.survey.data <- map_dfr(c("Cluster Survey V1.dta", "Cluster Survey 27.07.2016.dta"), ~read_cluster_survey_data(file.path(raw_data_path, .))) %>%
  group_by(cluster.id) %>%
  mutate(target.village.id = seq_len(n())) %>%
  ungroup %>%
  mutate(cluster.id = as.character(cluster.id),
         target.lat = ifelse(target.gps_work3 == 1 | target.gps_work3 == "Yes", target.gps3latitude, target.manual_lat3),
         target.lon = ifelse(target.gps_work3 == 1 | target.gps_work3 == "Yes", target.gps3longitude, target.manual_long3),
         found = factor(found, 1:4, labels = c("Found", "Name Changed", "Closed", "Not Found"))) %>%
  left_join(rct.cluster.selection %>%
              st_drop_geometry() %>%
              select(cluster.id, pot.lon, pot.lat),
            by = "cluster.id") 

cluster.survey.data %<>% filter(!is.na(cluster.id))




colnames(rct.schools.data)
rct.schools.data %>%
    select(contains("selected"))

# Create 1000m buffers around selected schools
anchor.buffers <- rct.schools.data %>%
  filter(!is.na(selected.targeted) & selected.targeted) %>%
  st_transform(kenya.proj4) %>%
  st_buffer(dist = 1000)

# Compute spatial containment for target villages
village_spatial_data <- cluster.survey.data %>%
  filter(!is.na(target.lon), !is.na(target.lat)) %>%
  st_as_sf(coords = c("target.lon", "target.lat"), crs = wgs.84, remove = FALSE) %>%
  st_transform(kenya.proj4) %>%
  mutate(
    # Find which buffer (if any) contains each village
    target.anchor.zone.cluster.id = {
      within_result <- st_within(., anchor.buffers, sparse = TRUE)
      sapply(within_result, function(x) {
        if (length(x) == 0) NA_character_ else anchor.buffers$pot.cluster.id[x[1]]
      })
    },
    target.in.anchor.zone = !is.na(target.anchor.zone.cluster.id) & cluster.id == target.anchor.zone.cluster.id
  ) %>%
  st_drop_geometry() %>%
  select(cluster.id, target.village.id, target.anchor.zone.cluster.id, target.in.anchor.zone)

# Join back to main data
cluster.survey.data %<>%
  left_join(village_spatial_data, by = c("cluster.id", "target.village.id")) %>%
  group_by(cluster.id) %>%
  mutate(num.valid.villages = sum(target.in.anchor.zone, na.rm = TRUE)) %>%
  ungroup

# Save output
# write_rds(cluster.survey.data, file.path(raw_data_path, "takeup_cluster_survey.rds"))

write_rds(
  cluster.survey.data, file.path(raw_data_path, "REATTEMPTED_takeup_cluster_survey.rds")
)