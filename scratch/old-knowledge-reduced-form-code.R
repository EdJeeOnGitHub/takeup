
disaggregate_beliefs = function(.data) {
  .data %>%
    mutate(
      doesnt_know_other_dewormed = obs_know_person - knows_other_dewormed,
      doesnt_think_other_knows = obs_know_person - thinks_other_knows
    ) %>%
    gather(variable, value,
      knows_other_dewormed_yes, knows_other_dewormed_no, doesnt_know_other_dewormed,
      thinks_other_knows_yes, thinks_other_knows_no, doesnt_think_other_knows
    ) %>%
    mutate(
      knowledge_type = case_when(
        str_detect(variable, "_yes") ~ "yes",
        str_detect(variable, "_no") ~ "no",
        str_detect(variable, "doesnt") ~ "doesn't know"
      ),
      belief_type = if_else(str_detect(variable, "think"), "2ord", "1ord"),
      prop = value / obs_know_person
    )
}

make_know_df = function(.data) {
  .data %>%
    filter(knowledge_type == "doesn't know") %>%
    mutate(prop_knows = 1 - prop) %>%
    group_by(cluster.id) %>%
    mutate(cluster_id = cur_group_id()) %>%
    ungroup() %>%
    mutate(
      county = factor(county),
      cluster.id = factor(cluster.id),
      assigned_treatment = assigned.treatment,
      signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
      signal = factor(signal, levels = c("no signal", "signal"))
    ) %>%
    left_join(
      cluster_expected_dist_df %>% mutate(cluster.id = as.numeric(cluster.id)),
      by = c("cluster_id" = "cluster.id")
    ) %>%
    mutate(
      standard_clust_expected_dist = clust_expected_dist / sd_of_dist,
      mu_d = standard_clust_expected_dist
    )
}
disagg_belief_all_df = all_data %>%
  mutate(female = gender == "female") %>%
  left_join(
    cov_analysis_data %>%
      select(cluster.id, mu_d, standard_clust_expected_dist) %>%
      mutate(cluster.id = as.numeric(cluster.id)) %>%
      unique(),
    by = "cluster.id"
  ) %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  select(any_of(belief_shared_vars), sms.treatment, dist.pot.group,
         all_of(l_cov_vars), standard_clust_expected_dist, contains("know")) %>%
  disaggregate_beliefs()

disagg_base_belief_data = cov_analysis_data %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  filter(obs_know_person > 0) %>%
  select(any_of(belief_shared_vars), dist.pot.group,
         all_of(l_cov_vars), standard_clust_expected_dist, contains("know")) %>%
  disaggregate_beliefs()

know_df = disagg_base_belief_data %>% make_know_df()
know_all_df = disagg_belief_all_df %>% make_know_df()


know_1_df = know_df  %>%
  filter(belief_type == "1ord") 
know_2_df = know_df  %>%
  filter(belief_type == "2ord")


know_1_all_df = know_all_df  %>%
  filter(belief_type == "1ord") %>%
  mutate(cluster_id = as.numeric(cluster.id)) 
