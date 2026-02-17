
create_balance_comparisons = function(fit) {
  comp_df = avg_comparisons(
    fit,
    variables = list("treat_dist" = "all")
  ) %>%
  as_tibble()


  comp_df = comp_df %>%
    mutate(
      lhs_treatment = str_extract(contrast, "treat: \\w+") %>% str_remove("treat: "),
      # remove first treat word and search after it for next treat word
      sub_str = str_extract(contrast, "(?<=treat).*"),
      rhs_treatment = str_extract(sub_str, "treat: \\w+") %>% str_remove("treat: "),

      lhs_dist = str_extract(contrast, "(?<=dist: )\\w+"),
      sub_str_dist = str_extract(contrast, "(?<=dist: ).*"),
      rhs_dist = str_extract(sub_str_dist, "(?<=dist: )\\w+")
    )  %>%
    select(-sub_str, -sub_str_dist)



  same_dist_subset_comp_df = comp_df  %>%
    filter(
      lhs_dist == rhs_dist,
      rhs_treatment == "control" | lhs_treatment == "control",
      lhs_treatment != rhs_treatment
    ) 

  same_dist_bra_cal_comp_df = comp_df %>%
    filter(lhs_dist == rhs_dist) %>%
    filter(str_detect(contrast, "bracelet") & str_detect(contrast, "calendar"))

    rhs_control_comp_df = same_dist_subset_comp_df %>%
      filter(rhs_treatment == "control") %>%
      bind_rows(
        same_dist_bra_cal_comp_df %>%
          filter(rhs_treatment == "calendar")
      )

    lhs_control_comp_df = same_dist_subset_comp_df %>%
      filter(rhs_treatment != "control") %>%
      bind_rows(
        same_dist_bra_cal_comp_df %>%
          filter(rhs_treatment == "bracelet")
      )

    lhs_control_comp_df = lhs_control_comp_df %>%
      mutate(
        new_estimate = estimate*-1, 
        new_statistic = statistic*-1, 
        new_conf.low = conf.high*-1,
        new_conf.high = conf.low*-1,
        new_lhs_treatment = rhs_treatment,
        new_rhs_treatment = lhs_treatment
      )  %>%
      mutate(
        estimate = new_estimate, 
        statistic = new_statistic,
        conf.low = new_conf.low,
        conf.high = new_conf.high, 
        lhs_treatment = new_lhs_treatment,
        rhs_treatment = new_rhs_treatment
      ) %>%
      select(-contains('new_'))

    rearranged_comp_df = bind_rows(
      lhs_control_comp_df, 
      rhs_control_comp_df
   ) %>%
   select(-contrast)


    control_mean_df = fit %>%
      tidy(conf.int = TRUE) %>%
      filter(str_detect(term, "control")) %>%
      mutate(
        lhs_treatment = "control", rhs_treatment = NA, 
        lhs_dist = if_else(str_detect(term, "close"), "close", "far"), 
        rhs_dist = lhs_dist
      ) %>%
      select(
        -term
      )
    sample_mean_df = tibble(
      lhs_treatment = "control",
      rhs_treatment = NA,
      lhs_dist = "combined",
      rhs_dist = "combined",
      estimate = fitstat(fit, type = "my", verbose = FALSE)$my
    )

  rearranged_comp_df = rearranged_comp_df %>%
    bind_rows(
      control_mean_df,
      sample_mean_df
    ) %>%
    mutate(comp_type = "treatment")

    ## Now within treatment across distances
    dist_control_mean_df = fit %>%
      tidy(conf.int = TRUE) %>%
      filter(str_detect(term, "close")) %>%
      mutate(
        lhs_treatment = str_extract(
          term, 
          "(?<=treat: )\\w+"),
        rhs_treatment = NA,
        lhs_dist = "close",
        rhs_dist = NA
      )

    dist_comp_df = comp_df %>%
      filter(
        rhs_dist != lhs_dist,
        lhs_treatment == rhs_treatment
      )  %>%
      select(-contrast) %>%
      bind_rows(
        dist_control_mean_df
      ) %>%
      mutate(comp_type = "distance")

    final_clean_comp_df = bind_rows(
      rearranged_comp_df,
      dist_comp_df
    )
    return(final_clean_comp_df)
}


perform_balance_joint_test = function(fit, var, joint_R, close_R, far_R) {
  county_0_mat = matrix(
    0,
    nrow = max(nrow(joint_R), nrow(close_R), nrow(far_R)),
    ncol = coef(fit) %>% length() - 8
    )

  resid_df = fixest::degrees_freedom(fit, type = "resid")
  close_test = car::lht(
    fit,
    cbind(close_R, county_0_mat[1:nrow(close_R), ]),
    error.df = resid_df,
    test = "F"
  )

  far_test = car::lht(
    fit,
    cbind(far_R, county_0_mat[1:nrow(far_R), ]),
    error.df = resid_df,
    test = "F"
  )

  joint_test = car::lht(
    fit,
    cbind(joint_R, county_0_mat[1:nrow(joint_R),]),
    error.df = resid_df,
    test = "F"
  )


  pvals = lst(
    joint_pval = joint_test$`Pr(>F)`[2],
    far_pval = far_test$`Pr(>F)`[2],
    close_pval = close_test$`Pr(>F)`[2]
  ) 

  return(pvals)
}

perform_balance_joint_test = function(fit, var, joint_R, close_R, far_R) {
  county_0_mat = matrix(
    0,
    nrow = max(nrow(joint_R), nrow(close_R), nrow(far_R)),
    ncol = coef(fit) %>% length() - 8
    )

  resid_df = fixest::degrees_freedom(fit, type = "resid")
  close_test = car::lht(
    fit,
    cbind(close_R, county_0_mat[1:nrow(close_R), ]),
    error.df = resid_df,
    test = "F"
  )

  far_test = car::lht(
    fit,
    cbind(far_R, county_0_mat[1:nrow(far_R), ]),
    error.df = resid_df,
    test = "F"
  )

  joint_test = car::lht(
    fit,
    cbind(joint_R, county_0_mat[1:nrow(joint_R),]),
    error.df = resid_df,
    test = "F"
  )


  pvals = lst(
    joint_pval = joint_test$`Pr(>F)`[2],
    far_pval = far_test$`Pr(>F)`[2],
    close_pval = close_test$`Pr(>F)`[2]
  ) 

  return(pvals)
}
n_variables = 8
# matrix R for test 
hyp_matrix = cbind(
  matrix(-1, nrow = n_variables - 1, ncol = 1 ), 
  diag(x = 1, nrow = n_variables - 1, ncol = n_variables)[, 1:(n_variables - 1)]
)

zero_matrix = matrix(0, nrow = 3, ncol = n_variables - 1) 
part_hyp_matrix = zero_matrix
for (i in 1:3) {
  part_hyp_matrix[i, 2*i] = 1
}

hyp_matrix_close = cbind(
  matrix(-1, nrow = 3, ncol = 1), 
  part_hyp_matrix
)

hyp_matrix_far = cbind(
  matrix(0, nrow = 3, ncol = 1),
  matrix(-1, nrow = 3, ncol = 1), 
  part_hyp_matrix[, 1:(ncol(part_hyp_matrix) - 1)]
)
