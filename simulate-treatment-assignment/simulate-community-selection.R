library(here)
library(tidyverse)
library(sf)
library(magrittr)

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"


# Identify valid target villages and clusters -----------------------------
rct.cluster.selection <- read_rds(here("data", "rct_cluster_selection_2.0.rds"))
rct.targetable.schools <- read_rds(here("data", "rct_targetable_schools_2.0.rds"))
rct_village_df = read_rds(here("data", "rct_target_villages_2.0.rds"))
# in rct_village_df, cluster.id and target.village.id tells us which were 
# actually selected from cluster.survey.data
rct_village_df = rct_village_df %>%
  mutate(
    ed_village_id = paste0(cluster.id, "_", target.village.id)
  )

rct_village_cw_df = rct_village_df %>%
  select(cluster.id, ed_village_id) %>%
  unique()


cluster.survey.data <- read_rds(here("data", "takeup_cluster_survey.rds"))
cluster.survey.data = cluster.survey.data %>%
  mutate(
    ed_village_id = paste0(cluster.id, "_", target.village.id)
  )

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
    set_colnames(cluster_villages$ed_village_id) %>%
    as_tibble() %>%
    mutate(
      dist.from.cluster.id = rct_cluster_selection$cluster.id,
      cluster.group = rct_cluster_selection$cluster.group
    )


  long_dist_mat = dist_mat %>%
    tidyr::gather(target.village.id, dist, -c(dist.from.cluster.id, cluster.group)) 
    
  sub_dist_mat = long_dist_mat %>%
    filter(dist <= 2500) %>%
    filter(dist != 0)

  return(sub_dist_mat)
}

rerandomize_village_targeting = function(seed, village_dist_mat, rct_cluster_selection) {
  set.seed(seed)
  rct_cluster_selection = rct_cluster_selection %>%
    as_tibble() %>%
    select(
      cluster.id
    ) %>%
    left_join(
      village_dist_mat,
      by = c("cluster.id" = "dist.from.cluster.id")
    )

  rct_cluster_selection = rct_cluster_selection %>%
    mutate(
      vill_dist_cat = case_when(
        dist <= 1250 ~ "close",
        dist > 1250 ~ "far"
      )
    ) %>%
    group_by(cluster.id) %>%
    mutate(
      random_dist_group = sample(c("close", "far"), 1),
      random_treat = case_when(
        cluster.group == "control.ink" ~ sample(c("control", "ink"), 1),
        cluster.group == "bracelet.airtime" ~ sample(c("bracelet", "calendar"), 1)
      ),
      no_dist_group_match = all(vill_dist_cat != random_dist_group)
    ) 

    match_treat_df = rct_cluster_selection %>%
      group_by(cluster.id) %>%
      filter(no_dist_group_match == FALSE) %>%
      filter(random_dist_group == vill_dist_cat) %>%
      sample_n(1)

    nonmatch_close_treat_df = rct_cluster_selection %>%
      filter(no_dist_group_match == TRUE) %>%
      filter(
        random_dist_group == "close"
      ) %>%
      group_by(cluster.id) %>%
      filter(
        dist == min(dist)
      )

      nonmatch_far_treat_df = rct_cluster_selection %>%
        filter(no_dist_group_match == TRUE) %>%
        filter(
          random_dist_group == "far"
        ) %>%
        group_by(cluster.id) %>%
        filter(
          dist == max(dist)
        )
    
    treat_df = bind_rows(
      match_treat_df,
      nonmatch_close_treat_df,
      nonmatch_far_treat_df
    ) %>%
    ungroup() %>%
    select(-vill_dist_cat)
    return(treat_df)
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


sim_int = 10
rct_pots = rct.cluster.selection %>%
  as_tibble() %>%
  pull(cluster.id)

sim_pots = simmed_cluster_selection[[sim_int]] %>%
  as_tibble() %>%
  pull(cluster.id)

setdiff(sim_pots, rct_pots)

sim_dist_mat = create_village_dist_mat(
    cluster.survey.data,
    simmed_cluster_selection[[sim_int]]
)

rerandomize_village_targeting(sim_int, sim_dist_mat, simmed_cluster_selection[[sim_int]])

sim_dist_mats = map(simmed_cluster_selection, ~ create_village_dist_mat(cluster.survey.data, .x))



pot_seeds = 1:length(simmed_cluster_selection)
vill_seeds = 1:100

seed_grid = expand.grid(pot_seed = pot_seeds, vill_seed =  vill_seeds) %>%
    as_tibble() 

library(furrr)
plan(multicore)
sim_df = future_map2_dfr(
    seed_grid$pot_seed,
    seed_grid$vill_seed,
    ~ rerandomize_village_targeting(.y, sim_dist_mats[[.x]], simmed_cluster_selection[[.x]]) %>%
        mutate(
            pot_seed = .x,
            vill_seed = .y
        ),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
)
sim_df = sim_df %>%
  rename(
    pot.cluster.id = cluster.id
  )


clean_sim_df = sim_df %>%
  left_join(
    rct_village_cw_df,
    by = c("target.village.id" = "ed_village_id")
  )

clean_sim_df %>%
    write_csv(
        here::here("data", "simulated-counterfactual-treatment-assignment.csv")
    )

cell_expected_dist_df = clean_sim_df %>%
    group_by(
        random_treat,
        random_dist_group
    ) %>%
    summarise(
        dist = mean(dist)
    )  %>%
    ungroup() 

cell_expected_dist_df %>%
    write_csv(
        here::here("data", "cell_expected_dist_df.csv")
    )


cluster_expected_dist = clean_sim_df %>%
    group_by(
        cluster.id
    ) %>%
    summarise(
        dist = mean(dist)
    )  


cluster_expected_dist %>%
    write_csv(
        here::here("data", "cluster_expected_dist.csv")
    )


summ_cluster_expected_dist = clean_sim_df %>%
  group_by(cluster.id, random_dist_group) %>%
  summarise(
    dist = mean(dist),
    pct_no_match = mean(no_dist_group_match)
  ) 


summ_cluster_expected_dist %>%
  write_csv(
    here::here("data", "summ_cluster_expected_dist.csv")
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
        random_dist_group = dist.pot.group,
        random_treat = assigned.treatment
    ) %>% 
    ungroup() %>%
    mutate(
        across(c(random_dist_group, random_treat), str_to_title)
    ) %>%
    mutate(
        random_dist_group = factor(random_dist_group, levels = c("Close", "Far")),
        random_treat = factor(random_treat, levels = c("Control", "Ink", "Calendar", "Bracelet"))
    )

summ_sim_df = clean_sim_df %>%
    group_by(
        vill_seed,
        pot_seed,
        random_treat,
        random_dist_group
    ) %>%
    summarise(
        dist = mean(dist)
    ) %>%
    ungroup() %>%
    mutate(
        across(c(random_dist_group, random_treat), str_to_title)
    ) %>%
    mutate(
        random_dist_group = factor(random_dist_group, levels = c("Close", "Far")),
        random_treat = factor(random_treat, levels = c("Control", "Ink", "Calendar", "Bracelet"))
    )

summ_sim_df %>%
    ggplot(aes(
        x = dist/1000,
        fill = random_dist_group
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

library(tidyverse)

cluster_expected_dist = read_csv(
        here::here("data", "cluster_expected_dist.csv")
    )


cluster_ids_in_rct =  clean_sim_df %>%
  filter(!is.na(cluster.id)) %>%
  pull(cluster.id) %>%
  unique()

clean_sim_df %>%
  group_by(cluster.id) %>%
  filter(
    cluster.id %in% cluster_ids_in_rct[1:5]
  ) %>%
  ggplot(aes(
    x = dist/1000
  )) +
  geom_histogram(bins = 20)  +
  facet_wrap(~cluster.id)

analysis_data %>%
  select(contains('dist')) %>%
  colnames()

clean_sim_df_treat = clean_sim_df %>%
  left_join(
    analysis_data %>%
      select(
        cluster.id,
        assigned_dist_group = dist.pot.group
      ) %>%
      mutate(cluster.id = as.character(cluster.id)) %>%
      unique(),
    by = "cluster.id"
  ) %>%
  mutate(
    assigned_dist_group = fct_na_value_to_level(assigned_dist_group, "Backup/Unassigned Cluster"),
    assigned_dist_group = fct_recode(assigned_dist_group, "Assigned Close" = "close", "Assigned Far" = "far")
  )

library(ggridges)

p_assignment_density = clean_sim_df_treat %>%
  group_by(
    assigned_dist_group
  ) %>%
  ggplot(aes(
    x = dist/1000,
    fill = assigned_dist_group,
    y = assigned_dist_group
  )) +
  geom_density_ridges(alpha = 0.5)  +
  theme_ridges() +
  theme(
    legend.position = "bottom"
  )  +
  labs(
    x = "Distance to PoT (km)",
    y = "Counterfactual Assignment Density",
    fill = "Experimental Assignment"
  )

ggsave(
  plot = p_assignment_density,
  filename = here::here("temp-plots", "counterfactual-assignment-density.pdf"),
  width = 8,
  height = 6
)
stop()



p_assignment_density_all = clean_sim_df_treat %>%
  filter(assigned_dist_group != "Backup/Unassigned Cluster") %>%
  group_by(
    assigned_dist_group
  ) %>%
  ggplot(aes(
    x = dist/1000,
    fill = assigned_dist_group,
    y = assigned_dist_group
  )) +
  geom_density_ridges(alpha = 0.5)  +
  theme_ridges() +
  theme(
    legend.position = "bottom"
  )  +
  labs(
    x = "Distance to PoT (km)",
    y = "Counterfactual Assignment Density",
    fill = "Experimental Assignment"
  )

ggsave(
  plot = p_assignment_density_all,
  filename = here::here("temp-plots", "counterfactual-assignment-density-close-far.pdf"),
  width = 8,
  height = 6
)

p_assignment_density_dens = clean_sim_df_treat %>%
  filter(assigned_dist_group != "Backup/Unassigned Cluster") %>%
  group_by(
    assigned_dist_group
  ) %>%
  ggplot(aes(
    x = dist/1000,
    fill = assigned_dist_group,
  )) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Distance to PoT (km)",
    y = "Density",
    fill = "Experimental Assignment"
  ) +
  theme(legend.position = "bottom")

ggsave(
  plot = p_assignment_density_dens,
  filename = here::here("temp-plots", "counterfactual-assignment-density-close-far-dens.pdf"),
  width = 8,
  height = 6
)


summ_test_df = clean_sim_df_treat %>%
  group_by(cluster.id, assigned_dist_group) %>%
  summarise(
    mean_dist = mean(dist)
  ) %>%
  ungroup()

summ_ks_test = ks.test(
  summ_test_df %>%
    filter(assigned_dist_group == "Assigned Far") %>%
    pull(mean_dist),
  summ_test_df %>%
    filter(assigned_dist_group == "Assigned Close") %>%
    pull(mean_dist)
)

broom::tidy(summ_ks_test) %>%
  mutate(
    interpretation = "p-value < 0.05 would suggest that the two vectors are not the same (i.e. drawn from different distributions)"
  ) %>%
  write_csv(
    here::here("data", "counterfactual-treatment-assignment-ks-test.csv")
  )
