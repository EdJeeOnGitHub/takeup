library(tidyverse)


# Original randomization of villages/clusters - `assigned.dist.cat` is
# cluster-level randomization by design, `pot.cluster.id` is the cluster ID used
# in the randomization. `cluster.dist.cat` is the close/far/mixed classification
# based on number of schools in the area
design_vill_clust_df = read_rds("data/rct_targetable_schools_2.0.rds")
# Potentially swapped post village survey - `village.dist.cat` is the cluster level 
# swapped assignment
swapped_clust_df = read_rds("data/rct_target_villages_2.0.rds")
# Processed in takeup_field_notebook.Rmd
processed_clust_df = read_rds("data/takeup_processed_cluster_strat.rds")
# Actual analysis data
analysis_df = read_csv("temp-data/analysis-data.csv", guess_max = 1e6)




clean_design_clust_df = design_vill_clust_df %>%
    group_by(pot.cluster.id) %>%
    arrange(pot.cluster.id) %>%
    mutate(pot.cluster.id = as.numeric(pot.cluster.id)) %>%
    select(cluster.dist.cat, assigned.dist.cat) %>%
    distinct() %>%
    arrange(pot.cluster.id)


clean_swapped_clust_df = swapped_clust_df %>%
    mutate(cluster.id = as.numeric(cluster.id)) %>%
    select(
        cluster.id, 
        village.dist.cat
        ) %>%
    arrange(cluster.id)

clean_processed_clust_df = processed_clust_df %>%
    select(cluster.id, dist.pot.group, assigned.treatment) %>%
    arrange(cluster.id)


clean_analysis_df = analysis_df %>%
    select(cluster.id, dist.pot.group, assigned.treatment) %>%
    distinct() %>%
    arrange(cluster.id)



all_clust_df = clean_design_clust_df %>%
    mutate(
        clust_classification_cat = cluster.dist.cat,
        clust_randomization_cat = assigned.dist.cat
    ) %>%
    select(pot.cluster.id, clust_classification_cat, clust_randomization_cat) %>%
    full_join(
        clean_swapped_clust_df %>% mutate(swap_cat = village.dist.cat) %>% 
            select(cluster.id, swap_cat),
        by = c("pot.cluster.id" = "cluster.id")
    ) %>%
    full_join(
        clean_processed_clust_df %>% mutate(processed_cat = dist.pot.group) %>%
            select(cluster.id, processed_cat),
        by = c("pot.cluster.id" = "cluster.id")
    ) %>%
    full_join(
        clean_analysis_df %>% mutate(analysis_cat = dist.pot.group) %>%
            select(cluster.id, analysis_cat),
        by = c("pot.cluster.id" = "cluster.id")
    )  %>%
    ungroup()


mismatched_clust_df = all_clust_df %>%
    filter(processed_cat != clust_randomization_cat)

write_csv(mismatched_clust_df, "docs/randomization-mismatches.csv")

mismatch_counts = mismatched_clust_df %>%
    count(clust_randomization_cat, processed_cat, name = "n_clusters") %>%
    mutate(direction = paste(clust_randomization_cat, "->", processed_cat))

cat("\nMismatch summary:\n")
print(mismatch_counts)
cat("\nTotal mismatched clusters:", nrow(mismatched_clust_df), "\n")
