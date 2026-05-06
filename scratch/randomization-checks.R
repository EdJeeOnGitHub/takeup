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


# Total clusters in the design
total_clusters = nrow(all_clust_df)
cat("\nTotal clusters in design:", total_clusters, "\n")


# Total realised Close/Far
realised_close_clusters = all_clust_df %>%
    filter(analysis_cat == "close") %>%
    nrow()


realised_far_clusters = all_clust_df %>%
    filter(analysis_cat == "far") %>%
    nrow()
cat("\nRealised Close clusters:", realised_close_clusters, "\n")
cat("Realised Far clusters:", realised_far_clusters, "\n")


# Feasibility drops (158 -> 144) ------------------------------------------
# Reasons sourced from:
#   scratch/distance/distance-recreation-attempt.Rmd (lines 49-53)
#   rct-design-fieldwork/takeup_rct_target_villages.R (lines 6-7)

drop_reasons = tribble(
  ~cluster_id, ~reason_group,                    ~reason,
  277,         "fieldwork",                      "Too close to another cluster (503)",
  491,         "fieldwork",                      "Problematic urban cluster",
  492,         "fieldwork",                      "Problematic urban cluster",
  1,           "fieldwork",                      "Village dispute about PoT",
  678,         "fieldwork",                      "Hostile community member",
  737,         "fieldwork",                      "Data fabrication and medication theft",
  99,          "cluster survey stage",            "Problems at cluster survey stage",
  691,         "cluster survey stage",            "Problems at cluster survey stage",
  892,         "cluster survey stage",            "Problems at cluster survey stage",
  1293,        "cluster survey stage",            "Problems at cluster survey stage",
  239,         "cluster survey stage",            "Problems at cluster survey stage",
  201,         "distance to other PoT",           "No longer valid: too close to another PoT",
  853,         "distance to other PoT",           "No longer valid: too close to another PoT",
  1402,        "distance to other PoT",           "No longer valid: too close to another PoT"
)

dropped_clust_df = all_clust_df %>%
  filter(is.na(analysis_cat)) %>%
  left_join(drop_reasons, by = c("pot.cluster.id" = "cluster_id"))

cat("\n=== Dropped clusters (158 -> 144) ===\n")
print(
  dropped_clust_df %>%
    select(cluster_id = pot.cluster.id, original_cat = clust_randomization_cat, reason_group, reason) %>%
    arrange(cluster_id)
)

cat("\nDropped by original category:\n")
print(dropped_clust_df %>% count(clust_randomization_cat, name = "n_dropped"))


# Count breakdown: design -> drops -> changes -> final --------------------

design_close = sum(all_clust_df$clust_randomization_cat == "close", na.rm = TRUE)
design_far   = sum(all_clust_df$clust_randomization_cat == "far",   na.rm = TRUE)

retained = all_clust_df %>% filter(!is.na(analysis_cat))
dropped  = all_clust_df %>% filter( is.na(analysis_cat))

dropped_close = sum(dropped$clust_randomization_cat == "close", na.rm = TRUE)
dropped_far   = sum(dropped$clust_randomization_cat == "far",   na.rm = TRUE)

changed = retained %>% filter(as.character(analysis_cat) != as.character(clust_randomization_cat))
close_to_far = sum(changed$clust_randomization_cat == "close")
far_to_close = sum(changed$clust_randomization_cat == "far")

final_close = sum(retained$analysis_cat == "close")
final_far   = sum(retained$analysis_cat == "far")

cat("\n=== Count breakdown ===\n")
cat(sprintf("%-25s  Close  Far\n", "Stage"))
cat(sprintf("%-25s  %5d  %3d\n", "Design",             design_close, design_far))
cat(sprintf("%-25s  %5d  %3d\n", "Dropped",            -dropped_close, -dropped_far))
cat(sprintf("%-25s  %5d  %3d\n", "After drops",        design_close - dropped_close, design_far - dropped_far))
cat(sprintf("%-25s  %5d  %3d\n", "Changed (->far)",    -close_to_far, +close_to_far))
cat(sprintf("%-25s  %5d  %3d\n", "Changed (->close)",  +far_to_close, -far_to_close))
cat(sprintf("%-25s  %5d  %3d\n", "Final analysis",     final_close, final_far))

