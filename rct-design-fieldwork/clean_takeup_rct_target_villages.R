if (interactive()) {
  library(tidyverse)
  library(magrittr)
  library(sf)
  library(ggspatial)
  raw_data_path = here::here("data")
}

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

# Distance computation helper functions -----------------------------------

#' Convert village coordinates to spatial points
#' @param villages Data frame with target.lon, target.lat
#' @return sf object with villages as points in Kenya UTM projection
create_village_sf <- function(villages) {
  st_as_sf(villages, coords = c("target.lon", "target.lat"), crs = wgs.84, remove = FALSE) %>%
    st_transform(kenya.proj4)
}

#' Calculate distance matrix from all PoTs to villages
#' @param pot_sf sf object of potential treatment clusters
#' @param villages_sf sf object of target villages
#' @param village_ids Village IDs for column naming
#' @return tibble with dist.from.cluster.id and one column per village (distances in meters)
calculate_distance_matrix <- function(pot_sf, villages_sf, village_ids) {
  st_distance(pot_sf, villages_sf) %>%
    units::drop_units() %>%
    set_colnames(village_ids) %>%
    magrittr::extract(, !duplicated(colnames(.)), drop = FALSE) %>%
    tibble::as_tibble() %>%
    mutate(dist.from.cluster.id = pot_sf$cluster.id,
           cluster.group = pot_sf$cluster.group)
}

#' Reshape distance matrix to long format and separate own-cluster vs other-cluster distances
#' @param dist_matrix Wide-format distance matrix from calculate_distance_matrix()
#' @param target_cluster_id The cluster.id these villages belong to
#' @return Long-format tibble with dist.to.pot (own cluster) and dist.to.other.pot (nearest other cluster per group)
reshape_distance_data <- function(dist_matrix, target_cluster_id) {
  dist_matrix %>%
    mutate(cluster.id = target_cluster_id) %>%
    tidyr::gather(target.village.id, dist, -c(dist.from.cluster.id, cluster.id, cluster.group)) %>%
    # Self-join to separate own-cluster distance (dist.to.pot) from other-cluster (dist.to.other.pot)
    left_join(.,
              select(., dist.from.cluster.id, target.village.id, dist),
              by = c("cluster.id" = "dist.from.cluster.id", "target.village.id"),
              suffix = c(".to.other.pot", ".to.pot")) %>%
    # Keep only other-cluster distances
    filter(dist.from.cluster.id != cluster.id) %>%
    select(-dist.from.cluster.id) %>%
    # For each village, keep only the nearest other-cluster PoT within each cluster.group
    group_by(target.village.id, cluster.group) %>%
    filter(row_number(dist.to.other.pot) == 1) %>%
    ungroup() %>%
    # Spread cluster.group (bracelet.airtime, control.ink) into columns
    spread(cluster.group, dist.to.other.pot)
}

#' Apply village validation rules based on distance thresholds
#' @param distance_data Output from reshape_distance_data() with bracelet.airtime and control.ink columns
#' @return Data with valid.target.village flag (TRUE if meets minimum distance to other clusters)
validate_target_villages <- function(distance_data) {
  distance_data %>%
    mutate(valid.target.village = bracelet.airtime >= 3860 & control.ink >= 3000,
           target.village.id = as.integer(target.village.id))
}

#' Compute distance metrics for one cluster's villages
#' @param cluster_villages Data frame for one cluster with target.village.id, target.lon, target.lat, cluster.id
#' @param pot_sf sf object of all potential treatment clusters
#' @return Tibble with distance metrics and validity flags
compute_cluster_village_distances <- function(cluster_villages, pot_sf) {
  villages_sf <- create_village_sf(cluster_villages)

  calculate_distance_matrix(pot_sf, villages_sf, cluster_villages$target.village.id) %>%
    reshape_distance_data(cluster_villages$cluster.id[1]) %>%
    validate_target_villages()
}

# Identify valid target villages and clusters -----------------------------
rct.cluster.selection <- read_rds(file.path(raw_data_path, "rct_cluster_selection_2.0.rds")) %>% st_as_sf()
rct.targetable.schools <- read_rds(file.path(raw_data_path, "rct_targetable_schools_2.0.rds"))
rct.schools.data <- read_rds(file.path(raw_data_path, "takeup_rct_schools.rds")) %>% st_as_sf()
cluster.survey.data = read_rds(file.path(raw_data_path, "REATTEMPTED_takeup_cluster_survey.rds"))

# Join randomized distance category info to cluster survey data
cluster.survey.data = cluster.survey.data %>%
  left_join(
    rct.targetable.schools %>%
      filter(selected.targeted) %>%
      select(village.dist.cat, pot.cluster.id),
      by = c("cluster.id" = "pot.cluster.id")
  )



clusters.to.drop <- c("99", "691", "892", "1293", "239") %>% # Problems with clusters at cluster survey stage; dropped from study
  c("201", "853", "1402") # These are are study cluster that are no longer valid based on distance to other PoT

rct.cluster.selection %<>% magrittr::extract(!.$cluster.id %in% clusters.to.drop, )

rct.cluster.selection %<>% left_join(rct.schools.data %>%
                                      st_transform(wgs.84) %>%
                                      mutate(lon = st_coordinates(.)[, "X"],
                                             lat = st_coordinates(.)[, "Y"]) %>%
                                      st_drop_geometry() %>%
                                      select(cluster.id, lon, lat)) %>%
  rename(pot.lon = lon, pot.lat = lat)


rct.cluster.selection = rct.cluster.selection %>%
  st_drop_geometry() %>%
  left_join(cluster.survey.data %>%
              distinct(cluster.id, found, num.valid.villages, gpslatitude, gpslongitude, gps2latitude, gps2longitude) %>%
              filter(!duplicated(cluster.id)),
            by = "cluster.id") %>%
  mutate(alt.pot.lat = gps2latitude,
         alt.pot.lon = gps2longitude,
         pot.lat = ifelse(!is.na(gpslatitude), #& found %in% c("Found", "Name Changed"),
                          gpslatitude,
                          ifelse(!is.na(gps2latitude), #& found %in% c("Closed", "Not Found"),
                                 alt.pot.lat,
                                 pot.lat)),
         pot.lon = ifelse(!is.na(gpslongitude), #& found %in% c("Found", "Name Changed"),
                          gpslongitude,
                          ifelse(!is.na(gps2longitude), #& found %in% c("Closed", "Not Found"),
                                 alt.pot.lon,
                                 pot.lon))) %>%
  st_as_sf(coords = c("pot.lon", "pot.lat"), crs = wgs.84, remove = FALSE) %>%
  st_transform(kenya.proj4)

# Compute village distances and validity -----------------------------------
# For each cluster, calculate:
# 1. Distance from village to its own PoT (dist.to.pot)
# 2. Distance to nearest other-cluster PoT per treatment group (bracelet.airtime, control.ink)
# 3. Validity flag (valid.target.village = TRUE if far enough from other clusters)

cluster.survey.data = cluster.survey.data %>%
  select(cluster.id, target.village.id, target.lat, target.lon,
          target.village_name, target.pop_households, village.dist.cat) %>%
  filter(!is.na(target.lat), !is.na(target.lon)) %>%
  arrange(cluster.id, target.village.id) %>%
  left_join(
    # Split by cluster and compute distance metrics for each cluster's villages
    filter(., !is.na(target.lon), !is.na(target.lat)) %>%
      split(.$cluster.id) %>%
      purrr::map_dfr(~compute_cluster_village_distances(.x, rct.cluster.selection)),
    by = c("cluster.id", "target.village.id")
  ) %>%
  left_join(
    # Add cluster distance category from targetable schools
    filter(rct.targetable.schools, selected.targeted) %>%
      select(pot.cluster.id, cluster.dist.cat),
    by = c("cluster.id" = "pot.cluster.id")
  )



cluster.survey.data = cluster.survey.data %>%
  group_by(cluster.id) %>% {
    if (first(.$village.dist.cat) == "far") arrange(., desc(dist.to.pot)) else arrange(., dist.to.pot)
  } %>%
  filter(!duplicated(target.village_name)) %>%
  ungroup %>%
  mutate(target.actual.village.dist.cat = ifelse(dist.to.pot <= (2500/2), "close", "far"),
         target.valid.dist.to.pot = target.actual.village.dist.cat == village.dist.cat,
         target.in.range = dist.to.pot <= 2500,
         dist.switchable = !target.valid.dist.to.pot & valid.target.village,
         valid.target.village = valid.target.village & target.in.range,
         target.village.pop.size = ifelse(target.pop_households <= median(target.pop_households, na.rm = TRUE), "small", "big")) %>%
  filter(!is.na(target.village.pop.size)) %>%
  group_by(cluster.id) %>%
  mutate(cluster.pop.size.strata = unique(target.village.pop.size) %>% { ifelse(length(.) == 2, "mixed", .) },
         actual.cluster.dist.cat = unique(target.actual.village.dist.cat) %>% { ifelse(length(.) == 2, "mixed", .) }) %>%
  ungroup()


any_switchable_clusters = cluster.survey.data %>%
  group_by(cluster.id) %>%
  summarize(any.dist.switchable = any(dist.switchable)) %>%
  filter(any.dist.switchable) %>%
  pull(cluster.id) 

# 583 obs - VILLAGES
cluster.survey.data 

# 150 - POTS
rct.cluster.selection

# Whether a given cluster has any valid target village that would force it to
# switch distance categories (e.g., "close" cluster finding only "far" villages)
summ_cluster_switchable_df = cluster.survey.data %>%
            mutate(backup.cluster = FALSE) %>% # cluster.id %in% rct.backup.clusters$cluster.id) %>%
            group_by(cluster.id, backup.cluster, cluster.pop.size.strata, actual.cluster.dist.cat) %>%
            summarize(any.buffer.valid.village = any(valid.target.village),
            # ANY VILLAGE SWITCHABLE => CLUSTER SWITCHABLE
                      dist.switchable = any(dist.switchable)) %>%
            ungroup()

rct.cluster.selection = rct.cluster.selection %>%
 left_join(summ_cluster_switchable_df,
          by = "cluster.id") %>%
  left_join(rct.targetable.schools %>%
              filter(selected.targeted) %>%
              select(pot.cluster.id, village.dist.cat, cluster.dist.cat),
            by = c("cluster.id" = "pot.cluster.id"))

# For plottings later
swapped_rct_cluster_df = rct.cluster.selection %>%
  as_tibble() %>%
  select(
    cluster.id, dist.switchable,  village.dist.cat, cluster.dist.cat,
    backup.cluster
  ) %>%
  mutate(
    original.village.dist.cat = village.dist.cat,
    new.village.dist.cat = ifelse(dist.switchable & !backup.cluster, ifelse(village.dist.cat == "close", "far", "close"), village.dist.cat)
  ) %>%
  select(-village.dist.cat) %>%
  filter(original.village.dist.cat != new.village.dist.cat)

swapped_rct_cluster_ids = swapped_rct_cluster_df %>% pull(cluster.id)


cluster.survey.data %>%
  filter(cluster.id %in% swapped_rct_cluster_ids[1:1]) %>%
  select(
    cluster.id, village.dist.cat, target.actual.village.dist.cat, dist.to.pot, valid.target.village, dist.switchable
  )


rct.cluster.selection = rct.cluster.selection %>%
  mutate(
    original.village.dist.cat = village.dist.cat,
    # ENTIRE POT SWITCHES CATEGORY IF ANY VILLAGE IN THAT POT FORCES A SWITCH (AND NOT A BACKUP CLUSTER)
    village.dist.cat = ifelse(dist.switchable & !backup.cluster, ifelse(village.dist.cat == "close", "far", "close"), village.dist.cat),
    any.buffer.valid.village = (dist.switchable & !backup.cluster) | any.buffer.valid.village)



# Display RCT clusters ----------------------------------------------------
rct.cluster.selection %>%
  st_drop_geometry() %>%
  filter(!any.buffer.valid.village) %>%
  select(cluster.id, cluster.group) %>%
  left_join(cluster.survey.data) %>%
  select(cluster.id, matches("valid"), bracelet.airtime, control.ink, dist.to.pot, village.dist.cat, cluster.group, dist.switchable)

# Village selection -------------------------------------------------------

rct.villages <- cluster.survey.data %>%
  filter(valid.target.village, cluster.id %in% rct.cluster.selection$cluster.id) %>%
  group_by(cluster.id) %>%
  sample_n(1) %>%
  ungroup

# write_rds(rct.villages, "rct_target_villages_2.0.rds")
orig_rct.villages <- read_rds(file.path(raw_data_path, "rct_target_villages_2.0.rds"))



unused.villages <- cluster.survey.data %>% anti_join(rct.villages, by = c("cluster.id", "target.village.id"))

# RCT target villages info table ------------------------------------------
write.rct.village.info <- function(.data, .file.name) {
  .data %>%
    select(cluster.id, target.village_name, target.lon, target.lat,
           matches("boda"), matches("target.*elder"), matches("target.*chv")) %>%
    left_join(rct.cluster.selection %>%
                st_drop_geometry() %>%
                select(cluster.id, school.name, county, constituency, ward, pot.lon, pot.lat),
              by = "cluster.id") %>%
    arrange(as.integer(cluster.id)) %>%
    write_csv(.file.name, na = "")

}

rct.villages %>%
  colnames()


rct.villages %>%
  select(target.village.id)

stop()
# rct villages get saved here

rct.villages %>% write.rct.village.info("REATTEMPTED-takeup_rct_villages.csv")

# Plot RCT clusters and villages ------------------------------------------
# 1. Data Preparation using sf objects
map_data_sf <- rct.cluster.selection %>%
  filter(any.buffer.valid.village) %>%
  st_transform(wgs.84) %>%
  st_drop_geometry() %>%
  transmute(lon = pot.lon, lat = pot.lat, type = "PoT", cluster.id) %>%
  bind_rows(rct.villages %>%
              transmute(lon = target.lon, lat = target.lat, type = "Target Village", cluster.id)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = wgs.84, remove = FALSE)

# 3. Create the buffers
pot_sf <- map_data_sf %>% filter(type == "PoT")

buffers_sf <- bind_rows(
  pot_sf %>% st_transform(kenya.proj4) %>% st_buffer(3000) %>% mutate(cluster.group = "control.ink"),
  pot_sf %>% st_transform(kenya.proj4) %>% st_buffer(4000) %>% mutate(cluster.group = "bracelet.calendar")
) %>% 
  st_transform(wgs.84)

# 4. Final Plotting with ggspatial
ggplot() +
  # Adds OpenStreetMap tiles as a background (Free/Open Source)
  annotation_map_tile(type = "osm", zoomin = 0, progress = "none") +
  
  # Plot Polygons (Buffers)
  geom_sf(data = buffers_sf, aes(fill = cluster.group), alpha = 0.25, color = "black", size = 0.1) +
  
  # Plot Points
  geom_sf(data = map_data_sf, aes(color = type), shape = 1, stroke = 1.5, alpha = 0.8) +
  
  # Formatting
  scale_color_manual("", values = c("red", "darkgreen")) +
  scale_fill_manual("Cluster Group", values = c("control.ink" = "blue", "bracelet.calendar" = "orange")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


# Intercluster distance ---------------------------------------------------
rct.villages %>%
  left_join(rct.cluster.selection %>% st_drop_geometry() %>% select(cluster.id, cluster.group)) %>%
  group_by(cluster.group) %>%
  do(dist.mat = st_as_sf(., coords = c("target.lon", "target.lat"), crs = wgs.84, remove = FALSE) %>%
       st_transform(kenya.proj4) %>%
       { st_distance(rct.cluster.selection, .) %>% units::drop_units() })


library(tidyverse)
library(sf)
library(ggspatial)

swapped_rct_cluster_df %>%
  filter(cluster.id == "43") 
plot_switch_clusters = "43"

# any_switchable_clusters[1:1]
# Filter the PoTs (Points of Treatment) for these clusters
switched_pots_sf <- rct.cluster.selection %>%
  filter(cluster.id %in% plot_switch_clusters) %>%
  st_as_sf(coords = c("pot.lon", "pot.lat"), crs = wgs.84) %>%
  st_transform(kenya.proj4) 

# Filter the Villages selected for these clusters
switched_villages_sf <- cluster.survey.data  %>%
  filter(cluster.id %in% plot_switch_clusters) %>%
  st_as_sf(coords = c("target.lon", "target.lat"), crs = wgs.84) %>%
  st_transform(kenya.proj4)

selected_village_sf = rct.villages %>%
  filter(cluster.id %in% plot_switch_clusters) %>%
  st_as_sf(coords = c("target.lon", "target.lat"), crs = wgs.84) %>%
  st_transform(kenya.proj4)

# 2. Create Distance Zones 
# "Close" is defined as < 1.25km (2500/2)
# "Far" is defined as > 1.25km but < 2.5km

zone_close <- switched_pots_sf %>% 
  st_buffer(1250) %>% 
  mutate(zone = "Close Zone (<1.25km)")

zone_valid <- switched_pots_sf %>% 
  st_buffer(2500) %>% 
  mutate(zone = "Max Range (2.5km)")

# Combine for plotting and transform back to WGS84 for ggspatial
plot_zones <- bind_rows(zone_valid, zone_close) %>% st_transform(wgs.84)
plot_pots <- switched_pots_sf %>% st_transform(wgs.84)
plot_villages <- switched_villages_sf %>% st_transform(wgs.84)
plot_selected_villages <- selected_village_sf %>% st_transform(wgs.84)

# Create connection lines for visual clarity

# 3. Plot the "Switched" Clusters -----------------------------------------
ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0, progress = "none") +
  
  # Draw the Zones
  geom_sf(data = plot_zones, aes(color = zone), fill = NA, size = 1, linetype = "dashed") +
  
  # Draw the PoT (Center)
  geom_sf(data = plot_pots, color = "black", shape = 3, size = 3, stroke = 2) +
  
  # Draw the Switching Villages
  geom_sf(data = plot_villages, aes(fill = "Switched Village"), shape = 21, size = 3, color = "black") +

  # Draw the selected village (if different from the switched village)
  geom_sf(data = plot_selected_villages, aes(fill = "Selected Village"), shape = 22, size = 4, color = "black") +

  
  # Formatting
  scale_color_manual("Distance Zones", values = c("Close Zone (<1.25km)" = "red", "Max Range (2.5km)" = "blue")) +
  scale_fill_manual("", values = c("Potential Villages" = "orange", "Selected Village" = "green")) +
  labs(title = "Visualizing Distance Switches",
       caption = "PoT marked by (+)") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Swapped Clusters Visualization ----------------------------------------
# Show all clusters that switched distance categories, with original designation

library(tidyverse)
library(sf)
library(ggspatial)

# 1. Get all switched clusters with their original designations
switched_clusters <- rct.cluster.selection %>%
  filter(dist.switchable) %>%
  select(cluster.id, original.village.dist.cat, village.dist.cat,
         pot.lon, pot.lat, school.name)

cat("\n=== Switched Clusters Summary ===\n")
switched_clusters %>%
  st_drop_geometry() %>%
  mutate(switch = paste(original.village.dist.cat, "→", village.dist.cat)) %>%
  count(switch) %>%
  print()

# 2. Get ALL villages for switched clusters (not just selected ones)
switched_villages_data <- cluster.survey.data %>%
  filter(cluster.id %in% switched_clusters$cluster.id) %>%
  left_join(
    switched_clusters %>% st_drop_geometry() %>%
      select(cluster.id, original.village.dist.cat, school.name),
    by = "cluster.id"
  ) %>%
  mutate(
    village_category = case_when(
      dist.to.pot <= 1250 ~ "Close (<1.25km)",
      dist.to.pot <= 2500 ~ "Far (1.25-2.5km)",
      TRUE ~ "Out of Range (>2.5km)"
    ),
    village_category = factor(village_category,
                              levels = c("Close (<1.25km)",
                                        "Far (1.25-2.5km)",
                                        "Out of Range (>2.5km)")),
    was_selected = target.village.id %in% rct.villages$target.village.id
  )

# 3. Create comprehensive multi-panel plot
plot_data_sf <- switched_villages_data %>%
  st_as_sf(coords = c("target.lon", "target.lat"), crs = wgs.84)

pots_sf <- switched_clusters %>%
  st_transform(wgs.84)

# Create distance rings for each PoT
rings_list <- list()
for(i in 1:nrow(switched_clusters)) {
  pot_row <- switched_clusters[i, ]
  pot_utm <- st_transform(pot_row, kenya.proj4)

  rings_list[[i]] <- bind_rows(
    pot_utm %>% st_buffer(1250) %>% mutate(ring = "Close threshold (1.25km)"),
    pot_utm %>% st_buffer(2500) %>% mutate(ring = "Max valid range (2.5km)")
  ) %>%
    st_transform(wgs.84) %>%
    mutate(
      cluster.id = pot_row$cluster.id,
      original.village.dist.cat = pot_row$original.village.dist.cat,
      facet_label = paste0("Cluster ", cluster.id, " (", school.name, ")\n",
                          "Original: ", original.village.dist.cat,
                          " → Switched to: ", pot_row$village.dist.cat)
    )
}
rings_sf <- bind_rows(rings_list)

# Add facet labels to other data
plot_data_sf <- plot_data_sf %>%
  left_join(
    rings_sf %>%
      st_drop_geometry() %>%
      select(cluster.id, facet_label) %>%
      distinct(),
    by = "cluster.id"
  )

pots_sf <- pots_sf %>%
  mutate(facet_label = paste0("Cluster ", cluster.id, " (", school.name, ")\n",
                             "Original: ", original.village.dist.cat,
                             " → Switched to: ", village.dist.cat))

# 4. Create the plot
ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0, progress = "none") +

  # Distance rings
  geom_sf(data = rings_sf, aes(linetype = ring),
          fill = NA, color = "black", size = 0.8) +

  # PoT centers
  geom_sf(data = pots_sf, shape = 3, size = 4, stroke = 2, color = "red") +

  # Villages colored by actual distance category
  geom_sf(data = plot_data_sf,
          aes(fill = village_category,
              size = was_selected,
              shape = was_selected),
          color = "black", stroke = 0.5, alpha = 0.8) +

  # Facet by cluster
  facet_wrap(~facet_label, scales = "free") +

  # Formatting
  scale_fill_manual("Village Distance\nfrom PoT",
                    values = c("Close (<1.25km)" = "lightgreen",
                              "Far (1.25-2.5km)" = "orange",
                              "Out of Range (>2.5km)" = "red")) +
  scale_size_manual("", values = c("TRUE" = 3, "FALSE" = 2),
                    labels = c("TRUE" = "Selected village", "FALSE" = "Not selected")) +
  scale_shape_manual("", values = c("TRUE" = 23, "FALSE" = 21),
                     labels = c("TRUE" = "Selected village", "FALSE" = "Not selected")) +
  scale_linetype_manual("Distance Rings",
                       values = c("Close threshold (1.25km)" = "dashed",
                                 "Max valid range (2.5km)" = "solid")) +
  labs(title = "Clusters That Switched Distance Categories",
       subtitle = "Shows why each cluster had to switch: Original designation vs actual available villages",
       caption = "PoT marked by red (+). Selected village shown as diamond.",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        strip.text = element_text(size = 9))

# 5. Print summary table
cat("\n=== Switched Clusters Detail ===\n")
switched_villages_data %>%
  group_by(cluster.id, school.name, original.village.dist.cat, village_category) %>%
  summarise(n_villages = n(),
            n_selected = sum(was_selected),
            .groups = "drop") %>%
  arrange(cluster.id, village_category) %>%
  print(n = 100)

# Single cluster example for detail
target_cluster_id <- rct.cluster.selection %>%
  filter(dist.switchable) %>%
  slice(2) %>%
  pull(cluster.id)

print(paste("Plotting Cluster ID:", target_cluster_id))

# 2. Get the PoT (Center)
pot_sf <- rct.cluster.selection %>%
  filter(cluster.id == target_cluster_id) %>%
  st_as_sf(coords = c("pot.lon", "pot.lat"), crs = wgs.84)

# 3. Get ALL Potential Villages
all_villages_sf <- cluster.survey.data %>%
  filter(cluster.id == target_cluster_id) %>%
  mutate(
    status = case_when(
      target.village.id %in% rct.villages$target.village.id ~ "Selected Village",
      valid.target.village ~ "Valid (But Not Selected)",
      TRUE ~ "Invalid / Skipped"
    )
  ) %>%
  st_as_sf(coords = c("target.lon", "target.lat"), crs = wgs.84)

# We replicate the single PoT coordinate to match the number of villages
village_coords <- st_coordinates(all_villages_sf)
pot_coords <- st_coordinates(pot_sf)

pot_utm <- st_transform(pot_sf, kenya.proj4) 

rings_sf <- bind_rows(
  pot_utm %>% st_buffer(1250) %>% mutate(label = "Close (1.25km)"),
  pot_utm %>% st_buffer(2500) %>% mutate(label = "Max (2.5km)")
) %>%
  st_transform(wgs.84)

# 6. Plot
ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0, progress = "none") +
  
  # Distance Rings
  geom_sf(data = rings_sf, aes(linetype = label), fill = NA, color = "black", size = 0.8) +
  
  
  # PoT Center
  geom_sf(data = pot_sf, shape = 3, size = 5, stroke = 2, color = "black") +
  
  # Village Points
  geom_sf(data = all_villages_sf, aes(fill = status, size = status), shape = 21, color = "black") +
  
  # Formatting
  scale_fill_manual("Village Status", values = c("Selected Village" = "blue", 
                                                 "Valid (But Not Selected)" = "green", 
                                                 "Invalid / Skipped" = "red")) +
  scale_size_manual("Village Status", values = c("Selected Village" = 4, 
                                                 "Valid (But Not Selected)" = 3, 
                                                 "Invalid / Skipped" = 2)) +
  labs(title = paste("Cluster:", target_cluster_id, "- Selection Logic"),
       subtitle = "Blue = Selected, Red = Invalid/Skipped. \nCheck if the Red dots forced the Blue dot into the wrong ring.",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")

# analysis_data = read_csv("temp-data/analysis-data.csv")

# joined_df = analysis_data %>%
#   # filter(cluster.id %in% swapped_rct_cluster_ids) %>%
#   select(cluster.id, dist.pot.group, cluster.dist.to.pot) %>%
#   unique() %>%
#   left_join(
#     cluster.survey.data %>%
#       # filter(cluster.id %in% swapped_rct_cluster_ids) %>%
#       select(cluster.id, village.dist.cat)  %>%
#       mutate(cluster.id = as.integer(cluster.id)) %>%
#       unique()
#   )  %>%
#   left_join(
#     rct.cluster.selection %>%
#       as_tibble() %>%
#       mutate(cluster.id = as.numeric(cluster.id)) %>%
#       select(cluster.id, rct_clust_dist.pot.group = village.dist.cat) 
#   )

# rct.cluster.selection %>%
#       as_tibble()  %>%
#       select(village.dist.cat)


# joined_df %>%
#   count(dist.pot.group)

# joined_df %>%
#   count(village.dist.cat)

# joined_df %>%
#   count(rct_clust_dist.pot.group)

# joined_df %>%
#   mutate(ed = cluster.dist.to.pot <= 1250) %>%
#   count(ed,dist.pot.group)