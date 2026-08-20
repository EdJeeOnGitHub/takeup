
#### Functions for RF bootstrap ------------------------------------------------
# Estimate and Predict to generate ATEs

# Run independent bootstrap draws across the CPU allocation. Each draw sets its
# own seed inside the estimation function, so fork scheduling does not affect
# the results. Forking also lets the workers share the large analysis frames by
# copy-on-write on Linux rather than serializing them for every draw.
bootstrap_map_dfr <- function(draws, fun, cores = NULL) {
  fun <- purrr::as_mapper(fun)
  if (is.null(cores)) {
    cores <- suppressWarnings(as.integer(Sys.getenv("TAKEUP_THREADS", "1")))
  }
  if (is.na(cores) || cores < 1L) cores <- 1L
  cores <- min(cores, length(draws))
  if (.Platform$OS.type == "windows") cores <- 1L

  if (cores == 1L) {
    return(purrr::map_dfr(draws, fun, .progress = TRUE))
  }
  message("Bootstrap: ", length(draws), " draws across ", cores, " fork workers")
  quiet_fun <- function(draw) {
    old_options <- options(fixest.notes = FALSE, dplyr.summarise.inform = FALSE)
    on.exit(options(old_options), add = TRUE)
    fun(draw)
  }
  results <- parallel::mclapply(
    draws, quiet_fun,
    mc.cores = cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("Parallel bootstrap worker failed: ", as.character(results[[which(failed)[1L]]]))
  }
  dplyr::bind_rows(results)
}

# Mark a linear discrete-distance specification as eligible for the cached
# weighted-cross-product bootstrap. The original function remains attached and
# is always used for the realised fit and as a validation/fallback path.
enable_fast_wls <- function(f, outcome, design_terms, required) {
  attr(f, "takeup_fast_wls") <- list(
    outcome = outcome,
    design_terms = design_terms,
    required = required
  )
  f
}

enable_fast_discrete_wls <- function(f, outcome, controls = character()) {
  enable_fast_wls(
    f, outcome,
    design_terms = c(
      "0 + assigned_treatment * assigned_dist_group", controls, "county"
    ),
    required = controls
  )
}

enable_fast_continuous_wls <- function(f, outcome, main_distance,
                                       interaction_distance = main_distance,
                                       controls = character()) {
  interactions <- sprintf(
    "I((assigned_treatment == '%s') * %s)",
    c("ink", "calendar", "bracelet"), interaction_distance
  )
  enable_fast_wls(
    f, outcome,
    design_terms = c(
      "0 + assigned_treatment", main_distance, interactions, controls, "county"
    ),
    required = unique(c(main_distance, interaction_distance, controls))
  )
}

prepare_fast_wls <- function(f, data, type = "APE") {
  enabled <- !identical(Sys.getenv("TAKEUP_FAST_WLS", "1"), "0")
  specification <- attr(f, "takeup_fast_wls", exact = TRUE)
  if (!enabled || is.null(specification) || type != "APE") return(NULL)

  required <- c(
    specification$outcome, "assigned_treatment", "assigned_dist_group",
    "county", "cluster.id", "cluster_id_rank", specification$required
  )
  if (!all(required %in% names(data))) return(NULL)

  prepared <- data
  prepared$assigned_treatment <- factor(
    as.character(prepared$assigned_treatment),
    levels = c("control", "ink", "calendar", "bracelet")
  )
  prepared$assigned_dist_group <- factor(
    as.character(prepared$assigned_dist_group), levels = c("close", "far")
  )
  prepared$county <- factor(prepared$county)

  rhs <- paste(specification$design_terms, collapse = " + ")
  design_formula <- as.formula(paste("~", rhs), env = environment(f))
  observed_matrix <- model.matrix(design_formula, prepared, na.action = na.pass)
  outcome <- prepared[[specification$outcome]]
  estimation_rows <- complete.cases(observed_matrix, outcome)
  if (qr(observed_matrix[estimation_rows, , drop = FALSE])$rank !=
      ncol(observed_matrix)) return(NULL)

  arms <- levels(prepared$assigned_treatment)
  prediction_data <- bind_rows(lapply(arms, function(arm) {
    mutate(prepared, assigned_treatment = factor(
      arm, levels = levels(prepared$assigned_treatment)
    ))
  }))
  prediction_matrix <- model.matrix(
    design_formula, prediction_data, na.action = na.pass
  )

  list(
    f = f,
    original_data = data,
    data = prepared,
    observed_matrix = observed_matrix,
    outcome = outcome,
    estimation_rows = estimation_rows,
    prediction_data = prediction_data,
    prediction_matrix = prediction_matrix
  )
}

fast_wls_bs_f <- function(seed, context, ...) {
  set.seed(seed)
  n_clusters <- length(unique(context$data$cluster.id))
  cluster_weights <- drop(generate_dirichlet(rep(1, n_clusters), 1))
  weights <- cluster_weights[context$data$cluster_id_rank]
  rows <- context$estimation_rows
  X <- context$observed_matrix[rows, , drop = FALSE]
  y <- context$outcome[rows]
  w <- weights[rows]

  coefficients <- tryCatch(
    solve(crossprod(X, w * X), crossprod(X, w * y)),
    error = function(e) NULL
  )
  if (is.null(coefficients)) {
    return(bayes_bs_f(seed, context$f, context$original_data, ...))
  }
  predictions <- drop(context$prediction_matrix %*% coefficients)
  prediction_data <- context$prediction_data
  prediction_data$pred <- predictions
  output <- prediction_data %>%
    select(
      assigned_dist_group, assigned_treatment,
      standard_cluster.dist.to.pot, pred, any_of("sms_treatment")
    ) %>%
    create_bs_preds(...)
  output$seed <- seed
  output
}

validate_fast_wls <- function(context) {
  if (is.null(context)) return(invisible(TRUE))
  direct <- fast_wls_bs_f(1L, context)
  reference <- bayes_bs_f(1L, context$f, context$original_data)
  keys <- intersect(
    c("assigned_treatment", "assigned_dist_group", "signal", "sms_treatment"),
    names(direct)
  )
  comparison <- full_join(
    direct %>% select(all_of(keys), direct = mean_pred),
    reference %>% select(all_of(keys), reference = mean_pred),
    by = keys
  ) %>%
    mutate(difference = direct - reference)
  difference <- max(abs(comparison$difference), na.rm = TRUE)
  if (!is.finite(difference) || difference > 1e-9) {
    print(comparison)
    stop("Fast WLS validation failed; maximum prediction difference = ",
         format(difference, scientific = TRUE))
  }
  message("Fast WLS validated; maximum prediction difference = ",
          format(difference, scientific = TRUE))
  invisible(TRUE)
}

bootstrap_regression_draws <- function(f, data, B_draws, type = "APE", ...) {
  context <- prepare_fast_wls(f, data, type)
  if (!is.null(context)) {
    validate_fast_wls(context)
    return(bootstrap_map_dfr(
      seq_len(B_draws),
      ~fast_wls_bs_f(seed = .x, context = context, ...)
    ))
  }
  bootstrap_function <- if (type == "APE") bayes_bs_f else bayes_bs_f_at_x
  bootstrap_map_dfr(
    seq_len(B_draws),
    ~bootstrap_function(seed = .x, f = f, data = data, ...)
  )
}

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
        data$wt = weights[data$cluster_id_rank]
    }
    if (any(is.na(data$wt))) {
      stop("NA weights found in pred_bs_f")
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
        data$wt = weights[data$cluster_id_rank]
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
    pvals = if_else(pvals == 0, "$<$0.001", as.character(pvals))
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
# Wrap a regression function f with Lee (2009) sample-trimming.
# base_data: full population (e.g. endline_data filtered to sms.control) with
#   columns assigned_treatment, assigned_dist_group, in_fob_sample, cluster.id.
# direction: "upper" drops lowest outcomes (maximises TE), "lower" drops highest.
# The outcome used for sorting must be named prop_knows in data.
make_lee_trim_f = function(f, base_data, direction = c("upper", "lower"), B_draws = NULL) {
  force(f)  # prevent lazy-eval recursion: without this, f in the closure resolves
            # to the reassigned (wrapped) f in create_regression_output, not the original
  direction = match.arg(direction)
  counter   = new.env(parent = emptyenv())
  counter$n = 0L
  function(data, ...) {
    counter$n = counter$n + 1L
    message(sprintf("Lee (%s) draw %d", direction, counter$n))
    # One weight per cluster — use first() to avoid floating-point distinct issues
    cluster_wt = data %>%
      group_by(cluster.id) %>%
      summarise(wt = first(wt), .groups = "drop")

    base_data_w = base_data %>%
      left_join(cluster_wt, by = "cluster.id", relationship = "many-to-one") %>%
      mutate(wt = replace_na(wt, 0))

    sel_rates = base_data_w %>%
      group_by(assigned_treatment, assigned_dist_group) %>%
      summarise(sel_rate = weighted.mean(in_fob_sample, wt, na.rm = TRUE), .groups = "drop")

    control_sel = sel_rates %>%
      filter(as.character(assigned_treatment) == "control") %>%
      select(assigned_dist_group, control_sel_rate = sel_rate)

    trim_fracs = sel_rates %>%
      left_join(control_sel, by = "assigned_dist_group") %>%
      mutate(
        assigned_treatment = as.character(assigned_treatment),
        assigned_dist_group = as.character(assigned_dist_group),
        trim_frac = pmax(0, (sel_rate - control_sel_rate) / sel_rate)
      ) %>%
      select(assigned_treatment, assigned_dist_group, trim_frac)


    data = data %>%
      select(-any_of("trim_frac")) %>%
      mutate(
        assigned_treatment  = as.character(assigned_treatment),
        assigned_dist_group = as.character(assigned_dist_group)
      ) %>%
      left_join(trim_fracs, by = c("assigned_treatment", "assigned_dist_group")) %>%
      group_by(assigned_treatment, assigned_dist_group) %>%
      arrange(if (direction == "upper") prop_knows else desc(prop_knows), .by_group = TRUE) %>%
      mutate(
        cum_wt_frac = cumsum(wt) / sum(wt),
        wt          = if_else(cum_wt_frac <= trim_frac, 0, wt)
      ) %>%
      ungroup()

    old_opts = options(fixest.notes = FALSE)
    on.exit(options(old_opts), add = TRUE)
    suppressWarnings(f(data, ...))
  }
}

create_regression_output = function(
                                    data, f,
                                    B_draws = getOption("takeup.bootstrap_draws", 500L),
                                    stat = params$stat,
                                    caption = "Average Treatment Effects: Reduced Form",
                                    dependent_var = "Dependent variable: Take-up",
                                    model_label = "Reduced Form",
                                    type = "APE",
                                    stars = TRUE,
                                    drop_H0s = FALSE,
                                    flip_calendar_sign = FALSE,
                                    lee_direction = NULL,
                                    lee_base_data = NULL
                                    ) {

  if (!is.null(lee_direction)) {
    f = make_lee_trim_f(f, lee_base_data, lee_direction, B_draws = B_draws)
  }

  actual_f <- if (type == "APE") {
    actual_bayesian_bs_fit
  } else {
    actual_bayesian_bs_fit_at_x
  }
  bs_draws <- bootstrap_regression_draws(f, data, B_draws, type = type)

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
    filter(assigned_treatment != "$|Calendar| - |Bracelet|$") %>%
    nice_kbl_table(
      cap = caption,
      outcome_var = dependent_var,
      model_label = model_label
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
      outcome_var = dependent_var,
      model_label = model_label
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

            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Any Signal = No Signal, $p$-value"  = "Signal Two-Side Pval"), 
            assigned_treatment = fct_recode(assigned_treatment, "$H0$: Bracelet = Calendar, $p$-value" = "Bracelet - Calendar Two-Side Pval")
        )  %>%
        filter(
          !(assigned_treatment %in% c("Bracelet - No Signal", "$H0$: Bracelet > Calendar, $p$-value", "$H0$: Any Signal > No Signal, $p$-value"))
        )
    
    return(tbl)
}

nice_kbl_table = function(tbl, cap, outcome_var = "Dependent variable: Take-up",
                          stat = params$stat, model_label = "Reduced Form") {
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
      setNames(4, model_label)
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
                            flip_calendar_sign = FALSE,
                            lee_direction = NULL, lee_base_data = NULL,
                            B_draws = getOption("takeup.bootstrap_draws", 500L)) {
  default_table_options = list(
    caption = "Average Treatment Effects: Reduced Form",
    dependent_var = "Dependent variable: Take-up",
    type = "APE",
    stars = TRUE,
    drop_H0s = FALSE,
    model_label = "Reduced Form"
  )
  table_options = modifyList(default_table_options, table_options)
  output = create_regression_output(
    data = data,
    f = regression_spec,
    caption = table_options$caption,
    dependent_var = table_options$dependent_var,
    model_label = table_options$model_label,
    type = table_options$type,
    stars = table_options$stars,
    drop_H0s = table_options$drop_H0s,
    stat = stat,
    flip_calendar_sign = flip_calendar_sign,
    lee_direction = lee_direction,
    lee_base_data = lee_base_data,
    B_draws = B_draws
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
