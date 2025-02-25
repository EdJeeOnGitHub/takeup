library(here)
library(tidyverse)
library(sf)
library(magrittr)

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"


# Identify valid target villages and clusters -----------------------------
rct.cluster.selection <- read_rds(here("data", "rct_cluster_selection_2.0.rds"))
rct.targetable.schools <- read_rds(here("data", "rct_targetable_schools_2.0.rds"))
cluster.survey.data <- read_rds(here("data", "takeup_cluster_survey.rds"))


simmed_cluster_selection = read_rds(here("data", "simmed_rct_clusters.rds"))


rct.schools.data = read_rds(here("data", "takeup_rct_schools.rds")) %>%
  st_as_sf(wgs.84) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2])



#' Create Village Dist Matrix
#'
#' Find the distance between each village in the cluster set with all the other
#' villages I think (?)
create_village_dist_mat = function(cluster_villages, rct_cluster_selection) {
  villages_sf = st_as_sf(cluster_villages,
    coords = c("target.lon", "target.lat"),
    crs = wgs.84
  ) %>%
    st_transform(kenya.proj4)
  dist_mat = st_distance(
    rct_cluster_selection %>% st_transform(kenya.proj4),
    villages_sf
  ) %>%
    units::drop_units()  %>%
    set_colnames(cluster_villages$target.village.id) %>%
    as_tibble() %>%
    mutate(
      dist.from.cluster.id = rct_cluster_selection$cluster.id,
      cluster.id = cluster_villages$cluster.id[1],
      cluster.group = rct_cluster_selection$cluster.group
    )
  long_dist_mat = dist_mat %>%
    tidyr::gather(target.village.id, dist, -c(dist.from.cluster.id, cluster.id, cluster.group))

  comb_long_dist_mat = left_join(
    long_dist_mat,
    long_dist_mat %>%
      select(dist.from.cluster.id, target.village.id, dist),
    by = c("cluster.id" = "dist.from.cluster.id", "target.village.id"),
    suffix = c(".to.other.pot", ".to.pot")
  ) %>%
    filter(dist.from.cluster.id != cluster.id) %>% 
    select(-dist.from.cluster.id) %>% 
    group_by(target.village.id, cluster.group) %>% 
    filter(row_number(dist.to.other.pot) == 1)  %>%
    ungroup() 

  wide_comb_dist_mat = comb_long_dist_mat %>%
    spread(cluster.group, dist.to.other.pot) %>%
    mutate(
      valid.target.village = bracelet.airtime >= 3860 & control.ink >= 3000,
      target.village.id = as.integer(target.village.id)
    )
  return(wide_comb_dist_mat)
}


randomize_village_targeting = function(seed, cluster_survey_data, rct_cluster_selection) {
  set.seed(seed)
  cluster_survey_data %>%
    filter(valid.target.village, cluster.id %in% rct_cluster_selection$cluster.id) %>%
    group_by(cluster.id) %>%
    sample_n(1) %>%
    ungroup() %>%
    select(
      cluster.id, target.village.id, dist.to.pot, village.dist.cat
    ) %>%
    left_join(
      rct_cluster_selection %>%
        as_tibble() %>%
        select(
          cluster.id,
          cluster.group
        ),
        by = "cluster.id"
    ) %>%
    mutate(
      random_treat = case_when(
        cluster.group == "control.ink" ~ sample(c("control", "ink"), 1),
        cluster.group == "bracelet.airtime" ~ sample(c("bracelet", "calendar"), 1)
      )
    )
}

randomize_village_targeting(2, cluster.survey.data, rct.cluster.selection)
randomize_village_targeting(2, cluster.survey.data, simmed_cluster_selection[[2]])

cluster_seeds = 1:100
pot_seeds = 1:length(simmed_cluster_selection)

seed_grid = expand.grid(cluster_seeds, pot_seeds) %>%
    as_tibble()

library(furrr)
plan(multicore)
sim_df = future_map2_dfr(
    seed_grid$Var1,
    seed_grid$Var2,
    ~ randomize_village_targeting(.x, cluster.survey.data, simmed_cluster_selection[[.y]]) %>%
        mutate(
            cluster_seed = .x,
            pot_seed = .y
        ),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
)

cell_expected_dist_df = sim_df %>%
    group_by(
        random_treat,
        village.dist.cat
    ) %>%
    summarise(
        dist = mean(dist.to.pot)
    )  %>%
    ungroup() 

cell_expected_dist_df %>%
    write_csv(
        here::here("data", "cell_expected_dist_df.csv")
    )


cluster_expected_dist = sim_df %>%
    group_by(
        cluster.id
    ) %>%
    summarise(
        dist = mean(dist.to.pot)
    )  

cluster_expected_dist %>%
    write_csv(
        here::here("data", "cluster_expected_dist.csv")
    )



load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)




monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


analysis_data <- monitored_nosms_data

summ_analysis_data = analysis_data %>%
    group_by(
        assigned.treatment,
        dist.pot.group
    ) %>%
    summarise(
        dist = mean(cluster.dist.to.pot)
    ) %>%
    rename(
        village.dist.cat = dist.pot.group,
        random_treat = assigned.treatment
    ) %>% 
    ungroup() %>%
    mutate(
        across(c(village.dist.cat, random_treat), str_to_title)
    ) %>%
    mutate(
        village.dist.cat = factor(village.dist.cat, levels = c("Close", "Far")),
        random_treat = factor(random_treat, levels = c("Control", "Ink", "Calendar", "Bracelet"))
    )

summ_sim_df = sim_df %>%
    group_by(
        cluster_seed,
        pot_seed,
        random_treat,
        village.dist.cat
    ) %>%
    summarise(
        dist = mean(dist.to.pot)
    ) %>%
    ungroup() %>%
    mutate(
        across(c(village.dist.cat, random_treat), str_to_title)
    ) %>%
    mutate(
        village.dist.cat = factor(village.dist.cat, levels = c("Close", "Far")),
        random_treat = factor(random_treat, levels = c("Control", "Ink", "Calendar", "Bracelet"))
    )

summ_sim_df %>%
    ggplot(aes(
        x = dist/1000,
        fill = village.dist.cat
    )) +
    facet_wrap(~random_treat) +
    geom_density(alpha = 0.5) +
    geom_vline(
        data = summ_analysis_data,
        aes(xintercept = dist/1000),
        linetype = "dashed"
    ) +
    theme_bw() +
    theme(
        legend.position = "bottom"
    ) +
    ggthemes::scale_fill_canva(
        "",
        palette = "Primary colors with a vibrant twist"
    ) +
    labs(
        x = "Distribution of Average Distance to PoT (km)",
        y = "Density",
        caption = "Average distance in realised experimental sample shown by dashed line."
    )
ggsave("temp-plots/simulated-dist-to-pot.pdf", width = 8, height = 6)
