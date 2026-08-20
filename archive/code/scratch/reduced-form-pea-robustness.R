
calculate_pea_at_points = function(data, close_point, far_point) {
  f = main_spec_regression
  f_signal = main_spec_signal_regression
  data$wt = 1
  fit = f(data) 
  pred_data = tibble(
    assigned_treatment = c("control", "bracelet", "calendar", "ink"),
  ) %>%
  expand(
    assigned_treatment,
    nesting(
      standard_cluster.dist.to.pot = c(close_point, far_point),
      assigned_dist_group = c("close", "far")
    ),
    county = "Kakamega"
  )

  pred = predict(fit, newdata = pred_data)
  pred
  pred_data$pred = pred
  pred_data = pred_data %>%
    group_by(assigned_dist_group) %>%
    mutate(
      te = pred - pred[assigned_treatment == "control"]
    )
  return(pred_data)
}


calculate_pea_at_points(analysis_data, 1.0, 3.0) 

max_dist = max(analysis_data$cluster.dist.to.pot)
min_dist = min(analysis_data$cluster.dist.to.pot)

close_points = seq(0, 1250, length.out = 25) / sd_of_dist
far_points = seq(1250, max_dist, length.out = 25) / sd_of_dist

dist_points = expand_grid(
  close_points = close_points,
  far_points = far_points
) %>%
  mutate(
    preds = map2(
      close_points,
      far_points,
      ~calculate_pea_at_points(analysis_data, .x, .y),
      .progress = TRUE
    )
  )




long_dist_points = dist_points %>%
  unnest(preds)  %>%
  mutate(
    dist = standard_cluster.dist.to.pot*sd_of_dist,
    close_points = round(close_points*sd_of_dist, 0),
    far_points = round(far_points*sd_of_dist, 0)
  )

te_df = 
long_dist_points %>%
  group_by(assigned_treatment, close_points, far_points) %>%
  mutate(
    te_diff = te[assigned_dist_group == "far"] - te[assigned_dist_group == "close"]
  ) %>%
  mutate(te_diff = if_else(assigned_dist_group == "close", NA_real_, te_diff)) %>%
  ungroup() %>%
  mutate(te = if_else(assigned_treatment == "control", NA_real_, te)) %>%
  select(
    assigned_treatment,
    dist,
    pred, te, te_diff,
    close_points,
    far_points
  )


te_df %>%
  pivot_longer(
    c(pred, te, te_diff),
    names_to = "type",
    values_to = "value"
  )  %>%
  filter(type != "te_diff") %>%
  mutate(
    type = case_when(
      type == "pred" ~ "Predicted Levels",
      type == "te" ~ "Treatment Effect",
    )
  ) %>%
  ggplot(aes(
    x = dist,
    y = value,
    colour = assigned_treatment
  )) +
  facet_wrap(~type) +
  geom_point() +
  labs(
    colour = "",
    x = "Distance to PoT (m)",
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


subset_close_points = close_points[seq(1, length(close_points), 4)]

dist_df = bind_rows(
analysis_data %>%
  group_by(
    assigned_treatment, assigned_dist_group
  ) %>%
  summarise(
    mean_dist = mean(cluster.dist.to.pot)
  ),
analysis_data %>%
  group_by(
     assigned_dist_group
  ) %>%
  summarise(
    mean_dist = mean(cluster.dist.to.pot)
  )
) 

wide_dist_df = dist_df %>%
  pivot_wider(
    names_from = assigned_dist_group,
    values_from = mean_dist
  )  %>%
  mutate(
    close_points = round(close, 0),
    far_points = round(far, 0),
    dist_diff = far - close,
    te_diff = 0
  ) %>%
  mutate(
    label = case_when(
      is.na(assigned_treatment) ~ str_glue(
        "Pooled"
      ),
      !is.na(assigned_treatment) ~ str_glue(
        "{assigned_treatment}"
      )
    )
  ) %>%
  mutate(
    standard_close = close/sd_of_dist,
    standard_far = far/sd_of_dist,
    standard_dist_diff = standard_far - standard_close
  )

wide_dist_df %>%
  mutate(
    preds = map2(
      standard_close,
      standard_far,
      ~calculate_pea_at_points(analysis_data, .x, .y),
      .progress = TRUE
    )
  ) %>%
  rename(orig_treatment = assigned_treatment) %>%
  unnest(preds) %>%
  filter(is.na(orig_treatment) | orig_treatment == assigned_treatment) %>%
  select(assigned_treatment, assigned_dist_group, standard_close, standard_far, pred) %>%
  mutate(dist_type = if_else(
    is.na(orig_treatment),
    "pooled",
    "cell"
  )) %>%
  select(dist_type, everything()) %>%
  ungroup()  %>%
  select(-orig_treatment) %>%
  mutate(dist_diff = standard_far - standard_close)  %>%
  group_by(assigned_dist_group, dist_type)  %>%
  mutate(
    te_treat_minus_control = pred - pred[assigned_treatment == "control"]
  ) %>%
  # filter(assigned_treatment == "bracelet" | assigned_treatment == "control") %>%
  group_by(assigned_treatment, dist_type) %>%
  mutate(
    level_far_minus_close = pred[assigned_dist_group == "far"] - pred[assigned_dist_group == "close"]
  ) %>%
  ungroup() %>%
  group_by(dist_type, assigned_treatment) %>%
  mutate(
    te_double_diff = te_treat_minus_control[assigned_dist_group == "far"] - te_treat_minus_control[assigned_dist_group == "close"]
  ) %>%
  ungroup() %>%
  arrange(dist_type, assigned_dist_group, assigned_treatment)  %>%
  select(
    -te_treat_minus_control,
    -level_far_minus_close,
    # -te_double_diff,
    -dist_diff
  ) %>%
  # arrange(assigned_treatment, assigned_dist_group, dist_type) %>%
  mutate(
    across(c(standard_close, standard_far), ~round(.x*sd_of_dist, 0))
  )
  



  select(-orig_treatment, -standard_close, -standard_far)  %>%

  pivot_wider(
    names_from = dist_type,
    values_from = pred
  )


te_df %>%
  filter(assigned_treatment != "control") %>%
  filter(close_points %in% round(subset_close_points*sd_of_dist, 0)) %>%
  mutate(dist_diff = far_points - close_points) %>%
  ggplot(aes(
    x = dist_diff,
    y = te_diff,
    colour = close_points
  )) +
  facet_wrap(~assigned_treatment, ncol = 1) +
  geom_point() +
  theme_bw() +
  labs(
    x = "Far - Close Distance Difference (m)",
    y = "Far - Close TE Difference",
    colour = "Close Position"
  )  



te_df %>%
  filter(assigned_treatment != "control") %>%
  filter(close_points %in% round(subset_close_points*sd_of_dist, 0))  %>%
  ggplot(aes(
    x = close_points,
    y = far_points,
    fill = te_diff
  )) +
  geom_tile() +
  facet_wrap(~assigned_treatment)  +
  geom_text(
    data = wide_dist_df %>%
      filter(assigned_treatment != "control"),
    aes(label = label),
    show.legend = FALSE,
    nudge_y = -40
  ) +
  geom_text(
    data = wide_dist_df %>%
      ungroup() %>%
      filter(is.na(assigned_treatment)) %>%
      select(-assigned_treatment),
    aes(label = label),
    show.legend = FALSE,
    nudge_y = 40
  )  +
  geom_point(
    data = wide_dist_df %>%
      filter(assigned_treatment != "control"),
      shape = 4
  ) +
  geom_point(
    data = wide_dist_df %>%
      ungroup() %>%
      filter(is.na(assigned_treatment)) %>%
      select(-assigned_treatment),
    shape = 5
  ) 

close_far_diff_df = long_dist_points %>%
  group_by(assigned_treatment, close_points, far_points) %>%
  summarise(
    dist_diff = unique(far_points - close_points),
    te_diff = te[assigned_dist_group == "far"] - te[assigned_dist_group == "close"]
  )


wide_dist_df

close_point_ends = c(close_points[1], close_points[length(close_points)])

close_far_diff_df %>%
  # filter(close_points %in% c(0, 1250)) %>%
  filter(assigned_treatment != "control") %>%
  ungroup() %>%
    filter(dist_diff >= 0) %>%
  ggplot(aes(
    x = dist_diff,
    y = te_diff,
    colour = assigned_treatment
  )) +
  geom_point() +
  facet_wrap(~close_points) +
  theme(legend.position = "bottom")

close_far_diff_df %>%
  filter(close_points %in% c(0, 1250)) %>%
  filter(assigned_treatment != "control") %>%
  ungroup() %>%
    filter(dist_diff >= 0) %>%
  ggplot(aes(
    x = dist_diff,
    y = te_diff,
    colour = assigned_treatment
  )) +
  geom_point() +
  facet_wrap(~close_points)



some_close_points = sample(close_points, 3)

close_far_diff_df %>%
  filter(close_points %in% some_close_points) %>%
  filter(assigned_treatment != "control") %>%
  ungroup() %>%
  mutate(
    close_points = round(close_points*sd_of_dist, 0),
    dist_diff = sd_of_dist*dist_diff/1000
    ) %>%
  ggplot(aes(
    x = dist_diff,
    y = te_diff,
    colour = assigned_treatment
  )) +
  geom_point() +
  facet_wrap(~close_points)




long_dist_points %>%
  mutate(
    dist_diff = far_points - close_points
  ) %>%
  # filter(close_points == 0) %>%
  ggplot(aes(
    x = dist_diff,
    y = te,
    colour = assigned_dist_group
  )) +
  geom_point() +
  facet_wrap(~assigned_treatment)

dist_points %>%
  unnest(preds) %>%
  ggplot(aes(
    x = standard_cluster.dist.to.pot,
    y = pred,
    colour = assigned_treatment
  )) +
  geom_point()

dist_points %>%
  unnest(preds) %>%
  filter(assigned_treatment == "bracelet") %>%
  ggplot(aes(
    x = close_points,
    y = far_points,
    fill = te
  )) +
  geom_tile() +
  facet_grid(assigned_dist_group~assigned_treatment)

ols_fit = feols(
  dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") | county,
  data = analysis_data,
  cluster = ~cluster.id
)


actual_fit = feglm(
  dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot + i(assigned_treatment, standard_cluster.dist.to.pot, "control") | county,
  family = "probit",
  data = analysis_data
)


analysis_data %>%
  group_by(
    assigned_treatment, assigned_dist_group
  ) %>%
  summarise(
    mean_dist = mean(cluster.dist.to.pot)
  ) 


long_dist_df = bind_rows(
analysis_data %>%
  group_by(
    assigned_treatment, assigned_dist_group
  ) %>%
  summarise(
    mean_dist = mean(cluster.dist.to.pot)
  )
  ,
analysis_data %>%
  group_by(
    assigned_dist_group
  ) %>%
  summarise(
    mean_dist = mean(cluster.dist.to.pot)
  )  %>%
  mutate(
    assigned_treatment = "pooled"
  )
) %>%
  mutate(
    standard_cluster.dist.to.pot = mean_dist/sd_of_dist
  )

pred_df = long_dist_df %>%
  filter(
    assigned_treatment == "pooled"
  ) %>%
  expand(
    assigned_treatment = c("control", "bracelet", "calendar", "ink"),
    standard_cluster.dist.to.pot
    )  %>%
    mutate(
      type = "pooled_ave"
    ) %>%
    bind_rows(
      long_dist_df %>%
        filter(
          assigned_treatment != "pooled"
        ) %>%
        select(
          assigned_treatment,  standard_cluster.dist.to.pot
        ) %>%
        mutate(type = "not_pooled")
    ) %>%
    ungroup()  %>%
    mutate(
      county = "Kakamega"
    )  %>%
  group_by(assigned_treatment, type) %>%
  mutate(
    assigned_dist_group = case_when(
      standard_cluster.dist.to.pot == min(standard_cluster.dist.to.pot) ~ "close",
      standard_cluster.dist.to.pot == max(standard_cluster.dist.to.pot) ~ "far"
    )
  ) %>%
  ungroup()

dist_pred_df = pred_df %>%
  expand(
    standard_cluster.dist.to.pot = seq(0.0, 4.0, 0.01),
    assigned_treatment,
    county
    ) 

dist_pred_df$pred  = dist_pred_df     %>%
  predict(
    actual_fit, newdata = .
  )


dist_pred_df %>%
  group_by(standard_cluster.dist.to.pot) %>%
  mutate(
    pred = pred - pred[assigned_treatment == "control"]
  ) %>%
  ggplot(aes(
    x = standard_cluster.dist.to.pot,
    y = pred,
    colour = assigned_treatment
  )) +
  geom_point(size = 3) +
  geom_vline(
    data = long_dist_df %>% 
    filter(assigned_treatment != "pooled"),
    aes(
      xintercept = standard_cluster.dist.to.pot,
      colour = assigned_treatment
    )
  ) +
  geom_vline(
    data = long_dist_df %>%
      filter(assigned_treatment == "pooled") %>%
      ungroup() %>%
        select(-assigned_treatment),
    aes(
      xintercept = standard_cluster.dist.to.pot
    ),
    colour = "black"
  ) +
  facet_wrap(~assigned_treatment)

dist_pred_df %>%
  ggplot(aes(
    x = standard_cluster.dist.to.pot,
    y = pred,
    colour = assigned_treatment
  )) +
  geom_point(size = 3) +
  geom_vline(
    data = long_dist_df,
    aes(
      xintercept = standard_cluster.dist.to.pot,
      colour = assigned_treatment
    )
  )

long_dist_df  


pred_df$pred = pred_df %>%
  predict(
    actual_fit, newdata = .
  )
pred_df %>%
  ggplot(aes(
    x = standard_cluster.dist.to.pot,
    y = pred,
    colour = assigned_treatment,
    shape = type
  )) +
  geom_point(size = 3)

pred_df %>%
  ggplot(aes(
    x = pred,
    y  = assigned_treatment,
    colour = type,
    shape = assigned_dist_group
  )) +
  geom_point(size = 3) +
  theme_bw()



  pivot_wider(
    names_from = assigned_dist_group,
    values_from = mean_dist
  ) %>%
  mutate(
    diff = far - close
  )