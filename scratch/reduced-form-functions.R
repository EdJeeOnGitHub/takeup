
#### Functions for RF bootstrap ------------------------------------------------
# Estimate and Predict to generate ATEs

# For a given set of IDs, create bs data - n.b. this allows cluster to appear 
# multiple times
create_bs_data = function(split_data, ids) {
    bs_df = do.call(rbind, split_data[ids])
    return(bs_df)
}


# Create bootstrapped predictions
create_bs_preds = function(pred_df, ...) {
    # check if combined already present for PEAs
    combined_present = "combined" %in% pred_df$assigned_dist_group
    preds = pred_df %>%
        group_by(
            assigned_treatment,
            assigned_dist_group,
            ...
        ) %>%
        summarise(
            mean_pred = mean(pred, na.rm = TRUE),
            .groups = "drop"
        )
    if (!combined_present) {
      preds = preds %>%
        bind_rows(
            pred_df %>%
                group_by(
                    assigned_treatment,
                    ...
                ) %>%
                summarise(
                    mean_pred = mean(pred, na.rm = TRUE),
                    .groups = "drop"
                ) %>%
                mutate(
                    assigned_dist_group = "combined"
                )
        )
    }
        

    signal_pred = pred_df %>%
        mutate(
          signal = case_when(
            assigned_treatment %in% c("bracelet", "ink") ~ "signal",
            assigned_treatment %in% c("calendar", "control") ~ "no signal"
          )
        ) %>%
        group_by(signal, assigned_dist_group, ...) %>%
        summarise(
            mean_pred = mean(pred, na.rm = TRUE),
            .groups = "drop"
        ) 
            
    if (!combined_present) {
      signal_pred = signal_pred %>%
            bind_rows(
            pred_df %>%
                mutate(
                  signal = case_when(
                    assigned_treatment %in% c("bracelet", "ink") ~ "signal",
                    assigned_treatment %in% c("calendar", "control") ~ "no signal"
                  )
                ) %>%
                group_by(
                    signal, 
                    ...
                ) %>%
                summarise(
                    mean_pred = mean(pred, na.rm = TRUE)
                ) %>%
                mutate(
                    assigned_dist_group = "combined"
                )
            )
    } 
          
    preds = bind_rows(preds, signal_pred)
       return(preds)
}



# Function to generate samples from a Dirichlet distribution
generate_dirichlet <- function(alpha, n) {
  # Generate matrix of gamma samples
  gamma_samples <- matrix(rgamma(n * length(alpha), shape = alpha, scale = 1), nrow = n)
  
  # Normalize rows to sum to 1
  dirichlet_samples <- gamma_samples / rowSums(gamma_samples)
  
  return(dirichlet_samples)
}

pred_bs_f = function(f, data, weights, realised_fit = FALSE) {
    if (realised_fit == TRUE) {
        data$wt = 1
    } else {
        data$wt = weights[data$cluster_id]
    }
    fit = f(data, weights = ~wt)

    pred_data = bind_rows(
      data %>%
        mutate(assigned_treatment = "bracelet"),
      data %>%
        mutate(assigned_treatment = "calendar"),
      data %>%
        mutate(assigned_treatment = "ink"),
      data %>%
        mutate(assigned_treatment = "control")
    )

    pred_data$pred = predict(fit, newdata = pred_data)
    data = pred_data %>%
        select(
            assigned_dist_group,
            assigned_treatment,
            standard_cluster.dist.to.pot,
            pred,
            any_of("sms_treatment")
        )


    return(data)
}

pred_bs_f_at_x = function(f,  data, weights, realised_fit = FALSE) {
    if (realised_fit == TRUE) {
        data$wt = 1
    } else {
        data$wt = weights[data$cluster_id]
    }
    fit = f(data, weights = ~wt)
    data = data %>%
        group_by(assigned_dist_group) %>%
        mutate(
          average_standard_cluster.dist.to.pot = mean(standard_cluster.dist.to.pot),
          average_dist.to.pot = mean(dist.to.pot)
        ) %>%
        ungroup() 


    collapsed_data = data %>%  
      mutate(
        standard_cluster.dist.to.pot = average_standard_cluster.dist.to.pot,
        dist.to.pot = average_dist.to.pot,
        county = "Kakamega"
      ) %>%
      select(
        standard_cluster.dist.to.pot,
        dist.to.pot,
        county,
        assigned_dist_group,
        assigned_treatment,
        signal,
        any_of("sms_treatment")
      ) %>%
      distinct() %>%
      bind_rows(
        data %>%
          group_by(assigned_treatment) %>%
          summarise(
              standard_cluster.dist.to.pot = mean(standard_cluster.dist.to.pot),
              dist.to.pot = mean(dist.to.pot),
              county = "Kakamega",
              signal = unique(signal)
          ) %>%
          mutate(
            assigned_dist_group = "combined"
          )
      )

    pea_pred = collapsed_data %>%  
      predict(
        fit, newdata = .
      )

    collapsed_data$pred = pea_pred


    pea_signal_pred = collapsed_data %>%
      predict(
        fit, newdata = .
      )

    collapsed_data$signal_pred = pea_signal_pred
    collapsed_data = collapsed_data %>%
        select(
            assigned_dist_group,
            assigned_treatment,
            signal,
            standard_cluster.dist.to.pot,
            pred,
            signal_pred,
            any_of("sms_treatment")
        )

    return(collapsed_data)

}

bayes_bs_f_at_x = function(seed, f, data, ...) {
    set.seed(seed)
    n_clusters = length(unique(data$cluster.id))
    alpha = rep(1, n_clusters)
    weights = generate_dirichlet(alpha, 1)
    bs_fit = pred_bs_f_at_x(f,  data, weights = weights) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 

bayes_bs_f = function(seed, f, data, ...) {
    set.seed(seed)
    n_clusters = length(unique(data$cluster.id))
    alpha = rep(1, n_clusters)
    weights = generate_dirichlet(alpha, 1)
    bs_fit = pred_bs_f(f, data, weights = weights) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 

actual_bayesian_bs_fit_at_x = function(seed, f,  data, ...) {
    bs_fit = pred_bs_f_at_x(f,  data, 1, realised_fit = TRUE) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    return(bs_fit)
} 
actual_bayesian_bs_fit = function(seed, f,  data, ...) {
    bs_fit = pred_bs_f(f,  data, 1, realised_fit = TRUE) %>%
        create_bs_preds(., ...)
    bs_fit$seed = seed
    data$wt = 1
    fit = f(data, weights = ~wt)
    return(list(bs_fit = bs_fit, fit = fit))
} 


add_predictions = function(draws, ...) {
    bra_minus_cal = draws %>%
            filter(assigned_treatment == "bracelet" | assigned_treatment == "calendar") %>%
            group_by(seed, assigned_dist_group, ...) %>%
            summarise(
                mean_pred = mean_pred[assigned_treatment == "bracelet"] - mean_pred[assigned_treatment == "calendar"],
            ) %>%
            mutate(
                assigned_treatment = "bracelet - calendar"
            )
      
    abs_cal_minus_abs_bra = draws %>%
            filter(assigned_treatment == "bracelet" | assigned_treatment == "calendar") %>%
            group_by(seed, assigned_dist_group, ...) %>%
            summarise(
                mean_pred = abs(mean_pred[assigned_treatment == "calendar"]) - abs(mean_pred[assigned_treatment == "bracelet"]),
            ) %>%
            mutate(
                assigned_treatment = "abs(calendar) - abs(bracelet)"
            )

    
    bra_minus_control = draws %>%
            filter(
              assigned_treatment == "bracelet" 
              ) %>%
            group_by(seed, assigned_dist_group, ...) %>%
            summarise(
                mean_pred = mean_pred[assigned_treatment == "bracelet"],
            ) %>%
            mutate(
                assigned_treatment = "bracelet - control"
            )

    bra_minus_no_signal = bind_rows(
      bra_minus_cal,
      bra_minus_control
    ) %>%
    group_by(seed, assigned_dist_group, ...) %>%
    summarise(
      mean_pred = mean(mean_pred)
    ) %>%
    mutate(assigned_treatment = "bracelet - no signal") 

    draws = bind_rows(draws, bra_minus_cal, abs_cal_minus_abs_bra, bra_minus_no_signal)

    fc_draws = draws %>%
    group_by(
        seed,
        assigned_treatment, ...
    ) %>%
    summarise(
        mean_pred = mean_pred[assigned_dist_group == "far"] - mean_pred[assigned_dist_group == "close"],
        assigned_dist_group = "far - close",
        .groups = "drop"
    )

    draws = bind_rows(draws, fc_draws) %>%
    ungroup() %>%
    mutate(
        assigned_treatment = factor(
            assigned_treatment,
            levels = c("control", "abs(calendar) - abs(bracelet)", "bracelet - calendar", "bracelet - no signal", "ink", "calendar", "bracelet")
        ),
        assigned_dist_group = factor(
            assigned_dist_group,
            levels = c("combined", "close", "far", "far - close")
        )
    ) 

    return(draws)
}

add_signal_predictions = function(draws, ...) {
    fc_draws = draws %>%
    group_by(
        seed,
        signal,
        ...
    ) %>%
    summarise(
        mean_pred = mean_pred[assigned_dist_group == "far"] - mean_pred[assigned_dist_group == "close"],
        assigned_dist_group = "far - close"
    )
    draws = bind_rows(draws, fc_draws) %>%
        ungroup() %>%
        mutate(
            signal = factor(
                signal,
                levels = c("no signal", "signal")
            ),
            assigned_dist_group = factor(
                assigned_dist_group,
                levels = c("combined", "close", "far", "far - close")
            )
        ) 

    return(draws)
}

create_tes = function(draws, ...) {
    draws %>%
        group_by(seed, assigned_dist_group, ...) %>%
        mutate(
            mean_pred = if_else(assigned_treatment == "control", mean_pred, mean_pred - mean_pred[assigned_treatment == "control"])
        )
}

create_cal_flip_tes = function(draws, ...){
    draws %>%
        group_by(seed, assigned_dist_group, ...) %>%
        mutate(
            mean_pred = case_when(
              assigned_treatment == "control" ~ mean_pred,
              assigned_treatment == "calendar" ~ -1*(mean_pred - mean_pred[assigned_treatment == "control"]),
              TRUE ~ mean_pred - mean_pred[assigned_treatment == "control"]
            ),
            cal_flipped = if_else(assigned_treatment == "calendar", TRUE, FALSE)
        )
}

create_signal_tes = function(draws, ...) {
    draws %>%
        group_by(seed, assigned_dist_group, ...) %>%
        mutate(
            mean_pred = if_else(signal == "no signal", mean_pred, mean_pred - mean_pred[signal == "no signal"])
        )
}


clean_signal_draws = function(draws, ...) {
    draws %>%
        filter(!is.na(signal)) %>%
        select(-assigned_treatment) %>%
        create_signal_tes(., ...) %>%
        add_signal_predictions(., ...) %>%
        rename(
            assigned_treatment = signal,
            estimate = mean_pred
        )
}

clean_te_flip_cal_draws = function(draws, ...) {
  draws %>%
        filter(!is.na(assigned_treatment)) %>%
        select(-signal) %>%
        create_cal_flip_tes(., ...) %>%
        add_predictions(., ...)  %>%
        rename(estimate = mean_pred)
}

clean_te_draws = function(draws, ...) {
    draws %>%
        filter(!is.na(assigned_treatment)) %>%
        select(-signal) %>%
        create_tes(., ...) %>%
        add_predictions(., ...)  %>%
        rename(estimate = mean_pred)
}



round_pval = function(pvals, digits = 3) {
    pvals = round(pvals, digits)
    pvals = if_else(pvals == 0, "<0.001", as.character(pvals))
    return(pvals)
}

add_summ_stats = function(bs_draws, actual_fit, ci_width = 0.95) {
    clean_tes = bs_draws %>%
      group_by(
          assigned_treatment,
          assigned_dist_group,
      ) %>%
      summarise(
          std_error = sd(estimate),
          conf.low = quantile(estimate, (1 - ci_width)/2, na.rm = TRUE),
          conf.high = quantile(estimate, 1 - (1 - ci_width)/2, na.rm = TRUE)
      ) %>%
      left_join(
          actual_fit,
          by = c("assigned_dist_group", "assigned_treatment")
      ) %>%
      mutate(
          pval = 2*pnorm(-abs(realised_pred)/std_error),
          oneside_pval = pnorm(-realised_pred/std_error)
      ) %>%
      mutate(
          pval = round_pval(pval, 3),
          oneside_pval = round_pval(oneside_pval, 3)
      ) %>%
      select(
          assigned_treatment, 
          assigned_dist_group, 
          realised_pred, 
          std_error, 
          conf.low,
          conf.high,
          pval, 
          oneside_pval) %>%
      rename(estimate = realised_pred) 
      return(clean_tes)
}

# wrapper function for all of the above
create_regression_output = function(data, f,  B_draws = 500, 
                                    stat = params$stat,
                                    caption = "Average Treatment Effects: Reduced Form",
                                    dependent_var = "Dependent variable: Take-up",
                                    type = "APE",
                                    stars = TRUE,
                                    drop_H0s = FALSE,
                                    flip_calendar_sign = FALSE
                                    ) {

  if (type == "APE") {
    bs_f = bayes_bs_f
    actual_f = actual_bayesian_bs_fit
  } else {
    bs_f = bayes_bs_f_at_x
    actual_f = actual_bayesian_bs_fit_at_x
  }
  bs_draws = map_dfr(
    1:B_draws,
    ~bs_f(
      seed = .x,
      f = f,
      data = data
    ),
    .progress = TRUE
    )

  if (flip_calendar_sign) {
    clean_te_draws_df = bs_draws %>%
      clean_te_flip_cal_draws()
  } else {  
  clean_te_draws_df = bs_draws %>%
    clean_te_draws()
  }

  clean_signal_draws_df = bs_draws %>%
    clean_signal_draws()

  realised_fit_output = actual_f(
    seed = "realised fit",
    f = f,
    data = data
  )
  realised_fit = realised_fit_output$bs_fit
  n_obs = nobs(realised_fit_output$fit)



  signal_fit = realised_fit %>%
    clean_signal_draws() %>%
    rename(realised_pred = estimate) %>%
    select(realised_pred, assigned_dist_group, assigned_treatment)

  if (flip_calendar_sign) {
    te_fit = realised_fit %>%
      clean_te_flip_cal_draws() %>%
      rename(realised_pred = estimate) %>%
      select(realised_pred, assigned_dist_group, assigned_treatment)
  } else {
    te_fit = realised_fit %>%
      clean_te_draws() %>%
      rename(realised_pred = estimate) %>%
      select(realised_pred, assigned_dist_group, assigned_treatment)
  }

  signal_summ = add_summ_stats(clean_signal_draws_df, signal_fit)
  te_summ = add_summ_stats(clean_te_draws_df, te_fit)

  pval_only_terms = c("bracelet - calendar", "signal")

  n_obs_df = tibble(
    assigned_treatment = factor("Observations"),
    pval = as.character(prettyNum(n_obs, big.mark = ",")),
    assigned_dist_group = unique(te_summ$assigned_dist_group),
    show_pval_only = TRUE,
    n_obs_line = TRUE
  )

  overall_summ = bind_rows(
    signal_summ,
    te_summ
  ) %>%
    mutate(
      show_pval_only = assigned_treatment %in% pval_only_terms
    ) %>%
    filter(assigned_treatment != "no signal")  %>%
    mutate(n_obs_line = FALSE) %>%
    bind_rows(n_obs_df) 

  if (drop_H0s) {
    overall_summ = overall_summ %>%
      filter(show_pval_only == FALSE | n_obs_line == TRUE)
  }

  default_tbl = overall_summ %>%
    prep_tbl(stat = stat, stars = stars) %>%
    nice_kbl_table(
      cap = caption,
      outcome_var = dependent_var
      )


  different_order_tbl = overall_summ %>%
    prep_tbl(stat = stat, stars = stars) %>%
    filter(assigned_treatment != "Observations") %>%
    mutate(
      assigned_treatment = fct_relevel(
        assigned_treatment, 
        c(
          "Control", 
          "Bracelet - No Signal", 
          "$|Calendar| - |Bracelet|$",
          "$H0$: Any Signal > No Signal, $p$-value",
          "$H0$: Any Signal $=$ No Signal, $p$-value",
          "$H0$: Bracelet > Calendar, $p$-value",
          "$H0$: Bracelet $=$ Calendar, $p$-value",
          "Ink", "Bracelet", "Calendar"
        ))) %>% 
      arrange(assigned_treatment) %>%
    nice_kbl_table(
      cap = caption,
      outcome_var = dependent_var
    )

  return(list(
    tidy_summary = overall_summ,
    default_tbl = default_tbl,
    different_order_tbl = different_order_tbl
  ))
}



### Table/Kable functions ------------------------------------------------------
prep_tbl = function(tes, stat = "ci", stars = FALSE) {

    tbl_dist_levels = c(
        "combined",
        "close",
        "far",
        "far - close"
    )

    tbl_contrast_levels = c(
        "bracelet",
        "calendar",
        "ink",
        "control",
        "Observations",
        "abs(calendar) - abs(bracelet)",
        "bracelet - no signal",
        "signal",
        "bracelet - calendar",
        "signal two-side pval",
        "bracelet - calendar two-side pval"
    )



    if (stat == "ci") {
        tes = tes %>%
            mutate(
                val = paste0(
                    "(",
                    round(conf.low, 3),
                    ", ",
                    round(conf.high, 3),
                    ")"
                )
            )
    } else if (stat == "std.error"){
        tes = tes %>%
            mutate(val = paste0("{(", round(std_error, 3), ")}"))
    } else {
        tes = tes %>%
            mutate(
                val = paste0("{(", round_pval(pval, 3), ")}")
            )
    }




    tbl =  tes %>%
        select(assigned_treatment, assigned_dist_group, estimate, conf.low, conf.high, val, pval, oneside_pval, show_pval_only, n_obs_line)  %>%
        mutate(
          show_stars = 
            ((assigned_treatment %in% c("bracelet", "calendar", "ink")) | assigned_dist_group == "far - close") & stars == TRUE,
          stars = case_when(
            pval < 0.001 ~ "***",
            pval < 0.05 ~ "**",
            pval < 0.1 ~ "*",
            TRUE ~ ""
          )
        ) %>%
        mutate(across(where(is.numeric), ~round(.x, 3))) %>%
        mutate(
            estim_std = case_when(
              show_pval_only == TRUE ~ linebreak(paste0(pval), align = "c"),
              show_stars == TRUE ~ linebreak(paste0(estimate, stars, "\n", str_glue("{val}")), align = "c"),
              TRUE ~ linebreak(paste0(estimate,"\n", str_glue("{val}")), align = "c")
            )
        ) %>%
        bind_rows(
            filter(., show_pval_only == TRUE & n_obs_line == FALSE) %>%
              mutate(
                estim_std = linebreak(paste0(pval), align = "c")
              ) %>%
              mutate(estim_std = as.character(estim_std)) %>%
              mutate(
                assigned_treatment = paste0(assigned_treatment, " two-side pval")
              )
        ) %>%
        select(assigned_treatment, assigned_dist_group, estim_std) %>%
        mutate(
            assigned_dist_group = factor(assigned_dist_group, tbl_dist_levels),
            assigned_dist_group = fct_relabel(assigned_dist_group, str_to_title),
            assigned_treatment = factor(assigned_treatment, tbl_contrast_levels)
        ) %>%
        arrange(assigned_dist_group, assigned_treatment) %>%
        pivot_wider(
            names_from = assigned_dist_group,
            values_from = estim_std
        ) %>%
        mutate(
            assigned_treatment = fct_relabel(assigned_treatment, str_to_title),
            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Any Signal > No Signal, $p$-value"  = "Signal"), 
            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Bracelet > Calendar, $p$-value" = "Bracelet - Calendar"),
            assigned_treatment = fct_recode(assigned_treatment, "$|Calendar| - |Bracelet|$" = "Abs(Calendar) - Abs(Bracelet)"),

            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Any Signal $\\neq$ No Signal, $p$-value"  = "Signal Two-Side Pval"), 
            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Bracelet $\\neq$ Calendar, $p$-value" = "Bracelet - Calendar Two-Side Pval")
        )  %>%
        filter(
          !(assigned_treatment %in% c("Bracelet - No Signal", "$H0$: Bracelet > Calendar, $p$-value", "$H0$: Any Signal > No Signal, $p$-value"))
        )
    
    return(tbl)
}

nice_kbl_table = function(tbl, cap, outcome_var = "Dependent variable: Take-up", stat = params$stat) {
  linesep_str = if_else(stat == "ci", "\\addlinespace", "")

  nice_kbl = tbl %>%
  kbl(
    col.names = c(
      # "Estimand", 
      # "Treatment", 
      outcome_var,
      paste0("(", 1:4, ")")
    ), 
    format = "latex", 
    linesep = linesep_str, 
    booktabs = TRUE, 
    escape = FALSE, 
    align = "lcccc", 
    caption = cap
  )  %>%
  kable_styling(
    latex_options = c("scale_down")
  ) %>%
  add_header_above(
    c(" ", 
      "Combined", 
      "Close", 
      "Far", 
      "Far - Close"
      ), 
    line = FALSE
  ) %>%
  add_header_above(
    c(
      " " = 1,
      "Reduced Form" = 4
      )
  ) %>%
  row_spec(c(3), hline_after = TRUE) 
}

custom_save_latex_table = function(table, table_name, table_output_path = params$table_output_path){
  table_conn = file(
    file.path(
      table_output_path, paste0(table_name, ".tex")
    )
  )

  attr(table, "kable_meta")$contents = str_replace_all(attr(table, "kable_meta")$contents, "removeme12345", " ")
  table[1] = str_replace_all(table[1], "removeme12345", " ")


  clean_table = table %>%
    str_remove(
      ., 
      fixed("\\begin{table}")
    ) %>%
    str_remove(
      .,
      "\\\\caption\\{.*\\}"
    ) %>%
    str_remove(
      ., 
      "\\\\end\\{table\\}"
    ) 
    
    
    clean_table %>%
      writeLines(
        table_conn
      )
    close(table_conn)

    return(table)
}

wrapper_function = function(data, regression_spec, tidy_summ_path, table_name, table_options = list(), stat = params$stat, 
                            flip_calendar_sign = FALSE) {
  default_table_options = list(
    caption = "Average Treatment Effects: Reduced Form",
    dependent_var = "Dependent variable: Take-up",
    type = "APE",
    stars = TRUE,
    drop_H0s = FALSE
  )
  table_options = modifyList(default_table_options, table_options)
  output = create_regression_output(
    data = data,
    f = regression_spec,
    caption = table_options$caption,
    dependent_var = table_options$dependent_var,
    type = table_options$type,
    stars = table_options$stars,
    drop_H0s = table_options$drop_H0s,
    stat = stat,
    flip_calendar_sign = flip_calendar_sign
  )
  output$tidy_summary %>%
    write_csv(tidy_summ_path)
  output$default_tbl %>%
    custom_save_latex_table(
      table_name = table_name
    )
  output$different_order_tbl %>%
    custom_save_latex_table(
      table_name = paste0(table_name, "_weird_order")
    )
    return(output)
}