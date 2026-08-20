#!/usr/bin/env Rscript
#
# Script: clean_hh_census_data.R
# Purpose: Generate cleaned household census data for TakeUp RCT
# Author: Generated from takeup_field_notebook.Rmd
# Date: 2026-02-10
#
# Description:
#   This script processes raw household census data from the TakeUp RCT,
#   performing the following operations:
#   - Loads and merges census data from multiple sources
#   - Fixes missing GPS coordinates
#   - Validates and corrects cluster assignments
#   - Identifies and removes invalid households
#   - Generates final cleaned dataset
#
# Usage:
#   Rscript clean_hh_census_data.R
#
# Output:
#   data/clean-data/clean-hh-census-data.csv
#
# Dependencies:
#   - tidyverse (dplyr, purrr, readr, stringr, etc.)
#   - sf (spatial features for geographic operations)
#   - here (path management)
#   - lubridate (date handling)
#   - glue (string interpolation)
#   - haven (for reading Stata .dta files)
#
# Note: This script uses modern sf package instead of deprecated rgeos/sp

# Load required packages --------------------------------------------------

library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(glue)
library(haven)

# Source helper functions (note: may still use sp/rgeos internally)
source(here("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))

# Constants ---------------------------------------------------------------

# Coordinate reference systems (using EPSG codes for sf)
# EPSG:4326 = WGS84 (lat/lon)
# EPSG:21036 = Kenya UTM Zone 36S
WGS_84_CRS <- 4326
KENYA_CRS <- 21036

# Date formats
DATETIME_FORMAT <- "%b %d, %Y %I:%M:%S %p"
TAKEUP_DATETIME_TYPE <- col_datetime(DATETIME_FORMAT)
TAKEUP_DATE_TYPE <- col_date(DATETIME_FORMAT)

# Clusters to exclude from analysis
CLUSTERS_TO_DROP <- c(
  277,  # Too close other cluster (I think 503)
  491, 492,  # Problematic urban clusters
  1,    # Village dispute about PoT
  678,  # Hostile community member
  737   # Data fabrication and medication theft
)

# UUIDs to manually drop from cluster corrections
MANUAL_CLUSTER_HH_TO_DROP <- c(
  "uuid:13ca71aa-d966-4956-8d8b-83673fd3f74a",
  "uuid:04d0d59b-e33b-49c6-b55f-ccf2eda14213",
  "uuid:fbedfc11-9e02-46c6-bf04-dfb817facefc",
  "uuid:611bf6c7-a84c-4210-809a-1d2b072dc682",
  "uuid:937590a7-b024-43c0-9bff-53d87defc663"
)

# Helper Functions --------------------------------------------------------

#' Get path to raw data file
raw_data_path <- function(filename) {
  here("data", "raw-data", filename)
}

#' Convert data frame to sf object
coords_to_sf <- function(data, crs = 4326) {
  data |>
    filter(!is.na(lon), !is.na(lat)) |>
    st_as_sf(coords = c("lon", "lat"), crs = crs, remove = FALSE)
}

#' Validate geographic coordinates
validate_coords <- function(data, bbox) {
  data |>
    mutate(
      invalid.coord =
        (!is.na(lon) & (lon > bbox["x", "max"] | lon < bbox["x", "min"])) |
        (!is.na(lat) & (lat > bbox["y", "max"] | lat < bbox["y", "min"]))
    )
}

#' Read and standardize TakeUp survey data
read_tu_data <- function(file_name,
                         submit_datetime_type = NULL,
                         other_types = NULL,
                         county_bbox) {
  col_types <- list(
    SubmissionDate = if (is.null(submit_datetime_type)) TAKEUP_DATETIME_TYPE else submit_datetime_type,
    manual_long = col_number(),
    manual_lat = col_number()
  ) |>
    c(other_types)

  read_csv(file_name, col_types = col_types, show_col_types = FALSE) |>
    mutate(isValidated = isValidated == "true") |>
    rename(
      lat = `gps-Latitude`,
      lon = `gps-Longitude`,
      cluster.id = cluster_id
    ) |>
    filter(deviceid != "(web)") |>
    mutate(
      lon = coalesce(lon, manual_long),
      lat = coalesce(lat, manual_lat)
    ) |>
    validate_coords(bbox = county_bbox)
}

#' Convert phone ownership codes to factor
factor_have_phone <- function(x) {
  factor(x, levels = 0:2, labels = c("No", "Yes", "Don't know number"))
}

#' Clean name columns
clean_names <- function(data) {
  data |>
    mutate(across(matches("1st|mid|2nd|clan"), str_to_upper)) |>
    mutate(across(
      matches("1st|mid|2nd|clan"),
      ~ if_else(.x %in% c("NONE", "N/A", "NA", "NO"), NA_character_, .x)
    )) |>
    mutate(across(
      matches("1st|2nd|clan"),
      ~ if_else(nchar(.x) < 2, NA_character_, .x)
    )) |>
    mutate(
      last_name = if_else(is.na(name2nd) & !is.na(name_mid), name_mid, name2nd),
      name_mid = if_else(is.na(name2nd) & !is.na(name_mid), NA_character_, name_mid)
    )
}

#' Identify closest cluster based on geographic proximity
identify_closest_cluster <- function(data, known_village_locations, kenya_crs = 21036) {
  # Convert household data to sf object
  hh_sf <- data |>
    filter(!is.na(lon), !is.na(lat)) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    st_transform(kenya_crs)

  # Convert village locations to sf object
  village_sf <- known_village_locations |>
    filter(!is.na(target.lon), !is.na(target.lat)) |>
    st_as_sf(coords = c("target.lon", "target.lat"), crs = 4326) |>
    st_transform(kenya_crs)

  # Calculate distance matrix
  dist_matrix <- st_distance(hh_sf, village_sf)

  # For each household, find closest village
  closest_info <- map_dfr(seq_len(nrow(hh_sf)), function(i) {
    distances <- as.numeric(dist_matrix[i, ])
    in_range <- distances <= 2000

    if (!any(in_range)) {
      tibble(
        KEY = data$KEY[i],
        min.dist = NA_real_,
        closest.cluster = NA_integer_
      )
    } else {
      min_dist <- min(distances[in_range])
      closest_idx <- which(distances == min_dist & in_range)[1]

      tibble(
        KEY = data$KEY[i],
        min.dist = min_dist,
        closest.village = village_sf$target.village.id[closest_idx],
        closest.cluster = village_sf$cluster.id[closest_idx]
      )
    }
  })

  closest_info |>
    select(KEY, min.dist, closest.cluster)
}

#' Get households outside boundary for a given cluster
get_outside_keys <- function(data, cluster_id, boundaries, kenya_crs = 21036) {
  # Get cluster boundary points
  boundary_points <- boundaries |>
    filter(cluster.id == cluster_id, boundary.type == "pre.census.convex.hull") |>
    filter(!is.na(lon), !is.na(lat))

  if (nrow(boundary_points) == 0) {
    return(character(0))
  }

  # Create convex hull and buffer
  cluster_boundary <- boundary_points |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    st_transform(kenya_crs) |>
    st_union() |>
    st_convex_hull() |>
    st_buffer(50) |>
    st_transform(4326)

  # Get households in this cluster
  hh_in_cluster <- data |>
    filter(cluster.id == cluster_id, !invalid.coord, !is.na(lon), !is.na(lat))

  if (nrow(hh_in_cluster) == 0) {
    return(character(0))
  }

  # Convert to sf and check which are outside boundary
  hh_sf <- hh_in_cluster |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326)

  # Find households NOT within boundary
  within_boundary <- st_within(hh_sf, cluster_boundary, sparse = FALSE)[, 1]
  outside_keys <- hh_in_cluster$KEY[!within_boundary]

  outside_keys
}

# Load reference data -----------------------------------------------------

message("Loading reference data...")

# County boundaries (for coordinate validation)
ke_lvl2_adm_data <- read_rds(here("data", "adm", "KEN_adm2.rds"))
rct_counties <- c("Busia", "Siaya", "Kakamega")
counties_adm_data <- ke_lvl2_adm_data[ke_lvl2_adm_data$NAME_1 %in% rct_counties, ]
county_bbox <- counties_adm_data@bbox

# Village codes
rct_village_codes <- read_csv(
  here("data", "village_codes_2.csv"),
  skip = 1,
  col_names = c("village.cluster.id", "village_name", "village"),
  show_col_types = FALSE
) |>
  mutate(village.cluster.id = as.integer(village.cluster.id))

# Cluster wave and county data
cluster_wave_county_data <- read_rds(here("data", "takeup_cluster_wave_county_5.0.rds"))

# All villages (for geographic matching)
rct_villages <- read_rds(here("data", "rct_target_villages_2.0.rds")) |>
  mutate(new.village = FALSE) |>
  bind_rows(
    read_rds(here("data", "rct_target_villages_2.0-4.rds")) |>
      mutate(new.village = TRUE)
  )

cluster_survey_data <- read_rds(here("data", "takeup_cluster_survey.rds"))

all_villages <- rct_villages |>
  bind_rows(anti_join(cluster_survey_data, rct_villages, by = c("cluster.id", "target.village.id"))) |>
  mutate(
    targeted.village = !is.na(new.village),
    cluster.id = as.integer(cluster.id),
    target.village_name = str_trim(target.village_name) |> str_replace_all("\\s+", " ")
  ) |>
  left_join(
    rct_village_codes,
    by = c("cluster.id" = "village.cluster.id", "target.village_name" = "village_name")
  ) |>
  left_join(cluster_wave_county_data, by = "cluster.id")

known_village_locations <- all_villages |>
  filter(targeted.village, !is.na(target.lon), !is.na(target.lat)) |>
  distinct(cluster.id, target.village.id, .keep_all = TRUE)

# Dictionaries for household and person IDs
hh_id_dict <- read_rds(here("data", "takeup_hh_id_dict.rds"))
person_id_dict <- read_rds(here("data", "takeup_person_id_dict.rds"))

# Pre-census data for boundary calculations
pre_census_data <- read_rds(here("data", "pre.census.processed.rds"))

# Get survey boundary for pre-census data
get_survey_boundary <- function(data, boundary_type, group_by = "cluster.id", kenya_crs = 21036) {
  data |>
    filter(!is.na(lon), !is.na(lat)) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    group_by(across(all_of(group_by))) |>
    summarise(geometry = st_union(geometry), .groups = "drop") |>
    st_convex_hull() |>
    st_transform(kenya_crs) |>
    st_buffer(dist = 0) |>  # Default buffer
    st_transform(4326) |>
    st_cast("POLYGON") |>
    mutate(
      boundary.type = boundary_type,
      group = paste(row_number(), boundary.type, sep = "-")
    )
}

pre_census_hhs_boundary <- pre_census_data |>
  mutate(cluster.id = as.integer(act.cluster.id)) |>
  get_survey_boundary("pre.census.convex.hull")

# Store as sf object for spatial operations
boundaries_sf <- pre_census_hhs_boundary

# Also create a data frame version for compatibility
boundaries <- pre_census_hhs_boundary |>
  st_coordinates() |>
  as_tibble() |>
  rename(lon = X, lat = Y) |>
  bind_cols(
    pre_census_hhs_boundary |>
      st_drop_geometry() |>
      slice(rep(seq_len(nrow(pre_census_hhs_boundary)),
                times = map_int(st_geometry(pre_census_hhs_boundary), ~ nrow(st_coordinates(.x)))))
  ) |>
  mutate(group = paste(group, boundary.type, sep = "-"))

# Clusters to drop from analysis
cluster_954_to_drop <- read_rds(here("data", "cluster_954_to_drop.rds"))

# Load Point-of-Treatment data --------------------------------------------

message("Loading Point-of-Treatment location data...")

# Cluster info entries to exclude
cluster_info_to_drop <- c(
  "uuid:3a491628-a0e8-4be0-9afc-c2cca08fc450",
  "uuid:fcd3137e-7ddb-4f2e-9140-680e7638e42c",
  "uuid:00fe16c6-5de1-4ff8-a9c8-92101ef4bd01",
  "uuid:e3a47c54-f7f2-4af0-bd5b-0de12d997cb0"
)

# Load cluster survey data (contains PoT locations)
cluster_info_v3 <- read_csv(
  raw_data_path("Cluster Survey V3 July 04.csv"),
  col_types = list(SubmissionDate = TAKEUP_DATE_TYPE),
  show_col_types = FALSE
)

cluster_info_v1 <- read_dta(here("data", "Cluster Survey V1.dta")) |>
  transmute(
    cluster.id = clusterid,
    SubmissionDate = parse_date(submissiondate, "%d/%m/%Y %T"),
    KEY = key,
    `gps2-Longitude` = gps2longitude,
    `gps2-Latitude` = gps2latitude,
    deviceid = as.character(deviceid),
    manual_long2, manual_lat2, location_type, alt_name, comments, enumerator
  )

cluster_info <- bind_rows(
  v3 = cluster_info_v3,
  v1 = cluster_info_v1,
  .id = "data.source"
) |>
  filter(
    !is.na(cluster.id),
    deviceid != "(web)",
    !KEY %in% cluster_info_to_drop
  ) |>
  arrange(cluster.id, SubmissionDate) |>
  group_by(cluster.id) |>
  slice_tail(n = 1) |>
  ungroup() |>
  rename(
    alt.pot.lon = `gps2-Longitude`,
    alt.pot.lat = `gps2-Latitude`
  ) |>
  mutate(
    alt.pot.lon = coalesce(alt.pot.lon, manual_long2),
    alt.pot.lat = coalesce(alt.pot.lat, manual_lat2),
    location_type = factor(
      location_type,
      levels = 1:5,
      labels = c("Clinic", "Church", "Market", "Home", "Other")
    ),
    cluster.id = as.integer(cluster.id)
  ) |>
  left_join(cluster_wave_county_data, by = "cluster.id")

# Load PoT verification data
pot_verify_data <- read_csv(
  raw_data_path("POT verification.csv"),
  col_types = list(SubmissionDate = TAKEUP_DATETIME_TYPE),
  show_col_types = FALSE
) |>
  rename(
    cluster.id = cluster_id,
    lon.verify = `gps-Longitude`,
    lat.verify = `gps-Latitude`
  ) |>
  select(-any_of("county")) |>  # Remove county column if it exists
  left_join(cluster_wave_county_data, by = "cluster.id") |>
  filter(
    wave == 1 | SubmissionDate >= "2016-10-16",
    !is.na(lon.verify), !is.na(lat.verify)
  ) |>
  group_by(cluster.id) |>
  mutate(num.entries = n()) |>
  slice_max(SubmissionDate, n = 1, with_ties = FALSE) |>
  ungroup()

# Merge cluster info with verification data
cluster_info %>%
  colnames()
pot_info <- cluster_info |>
  select(
    enumerator, enumerator_other, wave, county = county.x, cluster.id,
    alt.pot.lon, alt.pot.lat, alt_name, location_type,
    comments, SubmissionDate, data.source, KEY
  ) |>
  mutate(county = as.character(county)) |>
  left_join(
    pot_verify_data,
    by = c("wave", "county", "cluster.id"),
    suffix = c("", ".verify")
  )

# Calculate verification distance (original vs verified PoT coords)
pot_info_with_verify <- pot_info |>
  filter(
    !is.na(lon.verify), !is.na(lat.verify),
    !is.na(alt.pot.lon), !is.na(alt.pot.lat)
  ) |>
  mutate(
    verify.dist = map_dbl(seq_len(n()), ~ {
      pt1 <- st_point(c(alt.pot.lon[.x], alt.pot.lat[.x])) |>
        st_sfc(crs = WGS_84_CRS) |>
        st_transform(KENYA_CRS)
      pt2 <- st_point(c(lon.verify[.x], lat.verify[.x])) |>
        st_sfc(crs = WGS_84_CRS) |>
        st_transform(KENYA_CRS)
      as.numeric(st_distance(pt1, pt2)[1, 1])
    })
  )

pot_info <- pot_info_with_verify |>
  select(KEY, verify.dist) |>
  right_join(pot_info, by = "KEY")

# Create known PoT locations reference
known_pot_locations <- pot_info |>
  filter(
    !is.na(alt.pot.lon) | !is.na(lon.verify),
    !is.na(alt.pot.lat) | !is.na(lat.verify)
  ) |>
  mutate(
    alt.pot.lon = coalesce(alt.pot.lon, lon.verify),
    alt.pot.lat = coalesce(alt.pot.lat, lat.verify)
  ) |>
  distinct(cluster.id, .keep_all = TRUE)

# Load main census data ---------------------------------------------------

message("Loading main census data...")

# Pre-August 24th data
hh_census_data_pre24th <- read_tu_data(
  here("data", "census_gps_util23.csv"),
  submit_datetime_type = col_datetime("%d/%m/%Y %H:%M"),
  county_bbox = county_bbox
) |>
  filter(SubmissionDate >= "2016-08-22") |>
  mutate(across(c(present, return, number_individuals, instance, check), as.numeric))

# Main census data (post-August 24th)
hh_census_data <- read_tu_data(
  raw_data_path("Census.csv"),
  other_types = list(
    `gps-Latitude` = col_double(),
    `gps-Longitude` = col_double(),
    `gps-Altitude` = col_double(),
    `gps-Accuracy` = col_double(),
    deviceid = col_character(),
    gps_work = col_integer()
  ),
  county_bbox = county_bbox
) |>
  filter(SubmissionDate >= "2016-08-25") |>
  select(!starts_with("hhh_two_digits")) |>
  bind_rows(
    hh_census_data_pre24th |>
      select(!c(hhh_name1st_check, directions, starts_with("hhh_two_digits")))
  ) |>
  select(-c(RI, matches("^IC\\d{1,2}$"))) |>
  left_join(rct_village_codes, by = "village") |>
  mutate(
    village_name = factor(village_name),
    present = if_else(is.na(present) & SubmissionDate <= "2016-08-22", NA_real_, present),
    return = as.numeric(return),
    hhh_have_phone = factor_have_phone(hhh_have_phone)
  )

# Fix missing GPS data ----------------------------------------------------

message("Fixing missing GPS data...")

# Read various missing GPS files
missing_gps_data <- read_csv(
  here("data", "missing_gps.csv"),
  show_col_types = FALSE
) |>
  select(KEY, `gps-Latitude`, `gps-Longitude`, gps_work, manual_long, manual_lat) |>
  mutate(
    lon = coalesce(as.double(`gps-Longitude`), as.double(manual_long)),
    lat = coalesce(as.double(`gps-Latitude`), as.double(manual_lat)),
    data.source = "missing_gps"
  ) |>
  select(-starts_with("gps-"), -starts_with("manual_")) |>
  filter(!is.na(KEY), !is.na(lon), !is.na(lat))

missing_gps_149_data <- read_csv(
  here("data", "missing_gps_149.csv"),
  show_col_types = FALSE
) |>
  transmute(
    lon = gpslongitude,
    lat = gpslatitude,
    KEY = `key of form with missing GPS data`,
    data.source = "missing_gps_149"
  )

missing_gps_360_data <- read_csv(
  here("data", "missing_gps_360.csv"),
  show_col_types = FALSE
) |>
  transmute(
    lon = gpslongitude,
    lat = gpslatitude,
    KEY = `key of form missing GPS data`,
    data.source = "missing_gps_149"
  )

missing_gps_844_data <- read_csv(
  here("data", "missing_gps_844.csv"),
  show_col_types = FALSE
) |>
  rename(
    key.gps = `Key with GPS`,
    key.missing.gps = `Key missing GPS`
  ) |>
  inner_join(hh_census_data, by = c("key.gps" = "KEY")) |>
  transmute(lon, lat, KEY = key.missing.gps, data.source = "missing_gps_844")

# Clear GPS coords for households in missing_gps_844_data
hh_census_data <- hh_census_data |>
  mutate(
    mask = KEY %in% missing_gps_844_data$KEY,
    lon = if_else(mask, NA_real_, lon),
    lat = if_else(mask, NA_real_, lat)
  ) |>
  select(-mask)

# Combine all missing GPS data
missing_gps_data <- bind_rows(
  missing_gps_data,
  missing_gps_149_data,
  missing_gps_360_data,
  missing_gps_844_data
) |>
  validate_coords(bbox = county_bbox)

# Merge missing GPS data into main dataset
hh_census_data <- hh_census_data |>
  left_join(missing_gps_data, by = "KEY", suffix = c(".x", ".y")) |>
  filter(is.na(data.source) | data.source != "missing_gps_149" | cluster.id == 149) |>
  mutate(
    lon = coalesce(lon.x, lon.y),
    lat = coalesce(lat.x, lat.y),
    gps_work = coalesce(as.numeric(gps_work.y), gps_work.x),
    invalid.coord = coalesce(invalid.coord.y, invalid.coord.x)
  ) |>
  select(-matches("\\.[xy]$"), -data.source)

# Remove households marked for removal ------------------------------------

message("Removing households marked for exclusion...")

to_remove_files <- c(
  "takeup_returns_toremove.csv",
  "takeup_moved_notvalid_toremove.csv",
  "takeup_noconsent_toremove.csv"
)

hh_census_data <- hh_census_data |>
  anti_join(
    map_dfr(to_remove_files, ~ read_csv(here("data", .x), show_col_types = FALSE)),
    by = c("KEY" = "key")
  )

# Add cluster and household IDs -------------------------------------------

message("Adding cluster and household IDs...")

hh_census_data <- hh_census_data |>
  select(-county) |>
  left_join(cluster_wave_county_data, by = "cluster.id") |>
  left_join(hh_id_dict, by = "KEY")

# Process individual-level census data -----------------------------------

message("Processing individual-level census data...")

census_data <- read_csv(
  raw_data_path("Census-survey-individual.csv"),
  show_col_types = FALSE
) |>
  rename(KEY.individ = KEY) |>  # Rename before join
  inner_join(hh_census_data, by = c("PARENT_KEY" = "KEY")) |>
  rename(KEY = PARENT_KEY) |>  # Use PARENT_KEY as the household KEY
  mutate(
    two.digit.match = if_else(have_phone == 1, two_digits == two_digits_check, NA),
    have_phone = factor_have_phone(have_phone)
  ) |>
  group_by(KEY) |>
  mutate(
    num.individuals = n(),
    hh.has.phone = any(have_phone == "Yes"),
    hh.has.non.phone = any(have_phone != "Yes")
  ) |>
  ungroup() |>
  clean_names() |>
  filter(!is.na(name1st) & (!is.na(name_mid) | !is.na(name2nd))) |>
  left_join(person_id_dict, by = "KEY.individ")

# Merge household-level summaries back
hh_census_data <- hh_census_data |>
  inner_join(
    census_data |> distinct(KEY, num.individuals, hh.has.phone, hh.has.non.phone),
    by = "KEY"
  )

# Identify misidentified clusters -----------------------------------------

message("Identifying misidentified clusters...")

hh_census_data <- hh_census_data |>
  filter(
    !is.na(lon), !is.na(lat), !invalid.coord,
    cluster.id %in% known_village_locations$cluster.id
  ) |>
  identify_closest_cluster(known_village_locations, kenya_crs = 21036) |>
  right_join(hh_census_data, by = "KEY") |>
  mutate(misidentified.cluster = cluster.id != closest.cluster)

# Correct misidentified clusters ------------------------------------------

message("Correcting misidentified clusters...")

corrected_hhs <- hh_census_data |>
  filter(
    cluster.id %in% c(492, 254, 149, 360, 834, 587, 951, 1272, 735, 293, 1313),
    misidentified.cluster
  ) |>
  mutate(
    cluster.id = closest.cluster,
    misidentified.cluster = FALSE
  ) |>
  select(-c(wave, county)) |>
  left_join(cluster_wave_county_data, by = "cluster.id")

# Identify households to drop
cluster_hh_to_drop <- c(954, 197, 1442, 1242, 378, 431, 1426) |>
  map(~ get_outside_keys(hh_census_data, .x, boundaries, kenya_crs = 21036)) |>
  c(
    c(147, 587, 735, 293, 844) |>
      map(~ {
        hh_census_data |>
          filter(cluster.id == .x, is.na(closest.cluster)) |>
          pull(KEY)
      })
  ) |>
  unlist() |>
  c(MANUAL_CLUSTER_HH_TO_DROP)

# Apply corrections
hh_census_data <- corrected_hhs |>
  bind_rows(anti_join(hh_census_data, corrected_hhs, by = "KEY")) |>
  mutate(cluster.id = if_else(KEY %in% cluster_hh_to_drop, NA_integer_, cluster.id))

# Final data cleanup ------------------------------------------------------

message("Applying final data filters...")

hh_census_data <- hh_census_data |>
  filter(
    !invalid.coord,
    is.na(present) | (present == 1),
    is.na(moved) | (moved == 0) | (num.individuals > 0),
    is.na(consent) | (consent == 1),
    !is.na(cluster.id),
    !is.na(lon), !is.na(lat)
  ) |>
  mutate(
    old.cluster.id = if_else(cluster.id %in% CLUSTERS_TO_DROP, cluster.id, NA_integer_),
    cluster.id = if_else(cluster.id %in% CLUSTERS_TO_DROP, NA_integer_, cluster.id)
  ) |>
  group_by(cluster.id) |>
  filter(n() >= 50) |>
  ungroup()

# Calculate village-PoT distances -----------------------------------------

message("Calculating village centers and distances to PoT...")

# Calculate centroid of each cluster's households
village_centers <- hh_census_data |>
  filter(!is.na(cluster.id), !is.na(wave), !is.na(lon), !is.na(lat)) |>
  group_by(wave, county, cluster.id) |>
  summarise(
    # Calculate centroid using sf
    centroid_coords = list({
      st_as_sf(cur_data(), coords = c("lon", "lat"), crs = WGS_84_CRS) |>
        st_union() |>
        st_centroid() |>
        st_coordinates() |>
        as_tibble()
    }),
    .groups = "drop"
  ) |>
  unnest(centroid_coords) |>
  rename(lon = X, lat = Y) |>
  left_join(
    pot_info |>
      select(wave, county, cluster.id, starts_with("alt.pot."), ends_with(".verify")),
    by = c("wave", "county", "cluster.id")
  )

#' Calculate distances from village centers to PoTs
#' @param village_centers Data frame with village center coordinates
#' @param pot_coords_cols Character vector of length 2 with lon/lat column names
#' @param distance_threshold Distance in meters to categorize as "close" (default: 1250)
get_vill_pot_dist <- function(village_centers,
                               pot_coords_cols = c("alt.pot.lon", "alt.pot.lat"),
                               distance_threshold = 1250) {

  # Filter to rows with valid PoT coordinates
  valid_rows <- village_centers |>
    filter(
      !is.na(.data[[pot_coords_cols[1]]]),
      !is.na(.data[[pot_coords_cols[2]]])
    )

  if (nrow(valid_rows) == 0) {
    warning("No valid PoT coordinates found")
    return(village_centers |>
      mutate(
        dist.to.own.pot = NA_real_,
        dist.pot.group = NA_character_,
        closest.other = NA_real_
      ))
  }

  # Convert village centers to sf
  villages_sf <- valid_rows |>
    st_as_sf(coords = c("lon", "lat"), crs = WGS_84_CRS) |>
    st_transform(KENYA_CRS)

  # Convert PoT locations to sf using specified columns
  pots_sf <- valid_rows |>
    st_as_sf(coords = pot_coords_cols, crs = WGS_84_CRS, remove = FALSE) |>
    st_transform(KENYA_CRS)

  # Calculate distance matrix between all villages and all PoTs
  dist_matrix <- st_distance(villages_sf, pots_sf)

  # Extract diagonal (distance to own PoT)
  dist_to_own <- diag(dist_matrix)

  # Calculate distance to closest other PoT
  dist_matrix_no_diag <- dist_matrix
  diag(dist_matrix_no_diag) <- NA
  closest_other <- apply(dist_matrix_no_diag, 1, min, na.rm = TRUE)

  # Count number of PoTs within various distance thresholds
  distance_bins <- seq(distance_threshold, 7500, by = distance_threshold)
  count_within <- map(
    distance_bins,
    ~ rowSums(as.matrix(dist_matrix) < .x, na.rm = TRUE)
  ) |>
    set_names(paste0("num.within.", distance_bins))

  # Combine results
  valid_rows |>
    mutate(
      dist.to.own.pot = as.numeric(dist_to_own),
      dist.pot.group = factor(
        if_else(dist.to.own.pot <= distance_threshold, "close", "far"),
        levels = c("close", "far")
      ),
      closest.other = as.numeric(closest_other)
    ) |>
    bind_cols(as_tibble(count_within))
}

# Calculate distances using original PoT coordinates
vill_pot_dist <- village_centers |>
  filter(cluster.id %in% known_pot_locations$cluster.id) |>
  get_vill_pot_dist(pot_coords_cols = c("alt.pot.lon", "alt.pot.lat"))

# Calculate distances using verified PoT coordinates
verify_vill_pot_dist <- village_centers |>
  filter(
    cluster.id %in% known_pot_locations$cluster.id,
    !is.na(lon.verify), !is.na(lat.verify)
  ) |>
  get_vill_pot_dist(pot_coords_cols = c("lon.verify", "lat.verify"))

# Save output -------------------------------------------------------------

message("Saving cleaned household census data...")

# Create output directory if it doesn't exist
output_dir <- here("data", "clean-data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save main household census data
output_file <- file.path(output_dir, "clean-hh-census-data.csv")
write_csv(hh_census_data, output_file)

# Save village-PoT distance data as RData (for compatibility with original)
save(
  vill_pot_dist,
  verify_vill_pot_dist,
  village_centers,
  file = file.path(output_dir, "takeup_village_pot_dist.RData")
)


# Also save as individual CSV files for easier access
write_csv(vill_pot_dist, file.path(output_dir, "vill_pot_dist.csv"))
write_csv(verify_vill_pot_dist, file.path(output_dir, "verify_vill_pot_dist.csv"))
write_csv(village_centers, file.path(output_dir, "village_centers.csv"))

# Print summary -----------------------------------------------------------

message("\n=== Household Census Data Summary ===")
message(glue("Total households: {nrow(hh_census_data)}"))
message(glue("Total clusters: {n_distinct(hh_census_data$cluster.id, na.rm = TRUE)}"))
message(glue("Output saved to: {output_file}"))

# Print breakdown by wave
message("\nBreakdown by wave:")
hh_census_data |>
  count(wave, name = "n_households") |>
  mutate(wave = coalesce(as.character(wave), "Unknown")) |>
  arrange(wave) |>
  print()

message("\n=== Village-PoT Distance Summary ===")
message(glue("Village centers calculated: {nrow(village_centers)}"))
message(glue("Distances to original PoT: {nrow(vill_pot_dist)}"))
message(glue("Distances to verified PoT: {nrow(verify_vill_pot_dist)}"))

if (nrow(vill_pot_dist) > 0) {
  n_close <- sum(vill_pot_dist$dist.pot.group == "close", na.rm = TRUE)
  n_far <- sum(vill_pot_dist$dist.pot.group == "far", na.rm = TRUE)
  message(glue("Close villages (≤1250m): {n_close}"))
  message(glue("Far villages (>1250m): {n_far}"))

  message("\nDistance to own PoT summary (original coords):")
  print(summary(vill_pot_dist$dist.to.own.pot))
}

if (nrow(verify_vill_pot_dist) > 0) {
  message("\nDistance to own PoT summary (verified coords):")
  print(summary(verify_vill_pot_dist$dist.to.own.pot))
}

message("\n=== Output Files ===")
message(glue("  {output_file}"))
message(glue("  {file.path(output_dir, 'takeup_village_pot_dist.RData')}"))
message(glue("  {file.path(output_dir, 'vill_pot_dist.csv')}"))
message(glue("  {file.path(output_dir, 'verify_vill_pot_dist.csv')}"))
message(glue("  {file.path(output_dir, 'village_centers.csv')}"))

message("\n✓ Script completed successfully")
