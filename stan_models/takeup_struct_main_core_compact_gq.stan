// Compact generated quantities for takeup_struct_main_core.stan. The data and
// parameter blocks must remain schema-identical to the sampling model.

functions {
  #include struct_section_functions.stan
  // Heavy-tailed type robustness with bounded Newton solves; inactive in the
  // exact Gaussian branch.
  #include core_student_t_functions.stan
  #include core_finite_mixture_functions.stan
  #include core_asymmetric_observability_functions.stan

  // Probability and distance derivative for the selected recognition ladder.
  vector core_gq_recognition_channel(
      real distance,
      data int truth,
      data int treatment,
      data int recognition_structure,
      vector intercept,
      vector dist_slope,
      matrix arm_intercept,
      matrix arm_dist,
      data matrix contrast_basis) {
    vector[2] result = [1.0, 0.0]';
    if (recognition_structure < 2) {
      real eta = intercept[truth];
      real slope = 0;
      if (recognition_structure == 0) {
        eta += dot_product(arm_intercept[truth], contrast_basis[treatment]);
        slope = dist_slope[truth] +
          dot_product(arm_dist[truth], contrast_basis[treatment]);
        eta += slope * distance;
      }
      result[1] = inv_logit(eta);
      result[2] = result[1] * (1 - result[1]) * slope;
    }
    return result;
  }

  // First three entries are Pr(Yes, No, DK); last three are dPr/ddistance.
  vector core_gq_report_channel(
      real distance,
      data int truth,
      data int treatment,
      data int report_structure,
      data int is_public_signal,
      matrix report_intercept,
      matrix report_dist_slope,
      matrix report_arm_intercept,
      matrix report_arm_dist,
      data int report_arm_dist_hierarchical,
      matrix report_arm_dist_sd,
      vector definite_intercept,
      vector definite_dist_slope,
      matrix definite_arm_intercept,
      vector definite_public_signal_dist_slope,
      vector accuracy_intercept,
      matrix accuracy_arm_intercept,
      data matrix contrast_basis) {
    vector[6] result;
    if (report_structure == 2) {
      real definite_slope = definite_dist_slope[truth] +
        definite_public_signal_dist_slope[1] * is_public_signal;
      real definite_prob = inv_logit(
        definite_intercept[truth] +
        dot_product(definite_arm_intercept[truth], contrast_basis[treatment]) +
        definite_slope * distance
      );
      real accuracy_prob = inv_logit(
        accuracy_intercept[truth] +
        dot_product(accuracy_arm_intercept[truth], contrast_basis[treatment])
      );
      result[1:3] = core_two_stage_report_row(
        definite_prob, accuracy_prob, truth
      );
      result[4:6] = core_two_stage_report_row_derivative(
        definite_prob, definite_slope, accuracy_prob, truth
      );
    } else {
      vector[3] slope = rep_vector(0, 3);
      vector[2] logit;
      for (category in 1:2) {
        int first = (category - 1) * (cols(contrast_basis)) + 1;
        int last = category * cols(contrast_basis);
        real arm_slope = 0;
        if (report_structure == 0) {
          arm_slope = dot_product(
            report_arm_dist[truth, first:last], contrast_basis[treatment]
          );
          if (report_arm_dist_hierarchical) {
            arm_slope *= report_arm_dist_sd[truth, category];
          }
        }
        slope[category] = report_dist_slope[truth, category] + arm_slope;
        logit[category] = report_intercept[truth, category] +
          dot_product(
            report_arm_intercept[truth, first:last], contrast_basis[treatment]
          ) + slope[category] * distance;
      }
      result[1:3] = core_softmax_with_reference(logit[1], logit[2]);
      result[4:6] = result[1:3] .* (
        slope - rep_vector(dot_product(result[1:3], slope), 3)
      );
    }
    return result;
  }

  real core_gq_fixedpoint_residual(
      real cutoff,
      real benefit_cost,
      real mu_rep,
      real total_error_sd,
      real u_sd,
      data int use_u_in_delta) {
    array[1] int x_i = {use_u_in_delta};
    real delta;
    if (use_u_in_delta && u_sd > 0) {
      delta = expected_delta(
        cutoff, total_error_sd, u_sd, {0.0}, x_i
      );
    } else {
      delta = reputational_returns_normal(cutoff);
    }
    return cutoff + benefit_cost + mu_rep * delta;
  }

  // A bounded scalar bisection avoids sporadic algebra-solver failures in
  // extreme counterfactual draws. The paper's stable-equilibrium restriction
  // gives a unique crossing in this range for the retained posterior draws.
  real core_gq_find_fixedpoint_safe(
      real benefit_cost,
      real mu_rep,
      real total_error_sd,
      real u_sd,
      data int use_u_in_delta,
      data real alg_sol_rel_tol,
      data real alg_sol_f_tol,
      data real alg_sol_max_steps) {
    real lower_bound = -5.5;
    real upper_bound = 5.5;
    real lower_residual = core_gq_fixedpoint_residual(
      lower_bound, benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta
    );
    real upper_residual = core_gq_fixedpoint_residual(
      upper_bound, benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta
    );
    real midpoint = -benefit_cost;

    if (lower_residual * upper_residual <= 0) {
      for (step in 1:50) {
        real midpoint_residual;
        midpoint = 0.5 * (lower_bound + upper_bound);
        midpoint_residual = core_gq_fixedpoint_residual(
          midpoint, benefit_cost, mu_rep, total_error_sd, u_sd,
          use_u_in_delta
        );
        if (lower_residual * midpoint_residual <= 0) {
          upper_bound = midpoint;
          upper_residual = midpoint_residual;
        } else {
          lower_bound = midpoint;
          lower_residual = midpoint_residual;
        }
      }
    } else {
      // Defensive fallback for an unstable/no-crossing draw. Damped fixed-point
      // iteration stays finite, allowing downstream finite-draw checks to flag
      // it rather than losing the entire generated-quantities chain.
      for (step in 1:100) {
        real residual = core_gq_fixedpoint_residual(
          midpoint, benefit_cost, mu_rep, total_error_sd, u_sd,
          use_u_in_delta
        );
        midpoint = fmin(5.5, fmax(-5.5, midpoint - 0.5 * residual));
      }
    }
    return midpoint;
  }

  real core_gq_find_fixedpoint_student_t_safe(
      real benefit_cost,
      real mu_rep,
      real u_sd,
      data int use_u_in_delta,
      data real type_scale_sq,
      data vector type_precision,
      data vector type_weight) {
    real lower_bound = -8;
    real upper_bound = 8;
    real midpoint = -benefit_cost;
    vector[4] lower_moments = core_student_t_moments(
      lower_bound, u_sd, use_u_in_delta, type_scale_sq,
      type_precision, type_weight
    );
    vector[4] upper_moments = core_student_t_moments(
      upper_bound, u_sd, use_u_in_delta, type_scale_sq,
      type_precision, type_weight
    );
    real lower_residual = lower_bound + benefit_cost +
      mu_rep * lower_moments[3];
    real upper_residual = upper_bound + benefit_cost +
      mu_rep * upper_moments[3];
    if (lower_residual * upper_residual <= 0) {
      for (step in 1:60) {
        vector[4] midpoint_moments;
        real midpoint_residual;
        midpoint = 0.5 * (lower_bound + upper_bound);
        midpoint_moments = core_student_t_moments(
          midpoint, u_sd, use_u_in_delta, type_scale_sq,
          type_precision, type_weight
        );
        midpoint_residual = midpoint + benefit_cost +
          mu_rep * midpoint_moments[3];
        if (lower_residual * midpoint_residual <= 0) {
          upper_bound = midpoint;
          upper_residual = midpoint_residual;
        } else {
          lower_bound = midpoint;
          lower_residual = midpoint_residual;
        }
      }
    } else {
      for (step in 1:120) {
        vector[4] midpoint_moments = core_student_t_moments(
          midpoint, u_sd, use_u_in_delta, type_scale_sq,
          type_precision, type_weight
        );
        real residual = midpoint + benefit_cost +
          mu_rep * midpoint_moments[3];
        midpoint = fmin(8, fmax(-8, midpoint - 0.4 * residual));
      }
    }
    return midpoint;
  }

  real core_gq_find_fixedpoint_finite_mixture_safe(
      real benefit_cost, real mu_rep, real u_sd,
      data int use_u_in_delta, vector component_mean,
      vector component_variance, vector component_weight) {
    real lower_bound = -8;
    real upper_bound = 8;
    real midpoint = -benefit_cost;
    real lower_residual = lower_bound + benefit_cost + mu_rep *
      core_finite_mixture_moments(
        lower_bound, u_sd, use_u_in_delta, component_mean,
        component_variance, component_weight
      )[3];
    real upper_residual = upper_bound + benefit_cost + mu_rep *
      core_finite_mixture_moments(
        upper_bound, u_sd, use_u_in_delta, component_mean,
        component_variance, component_weight
      )[3];
    if (lower_residual * upper_residual <= 0) {
      for (step in 1:60) {
        real midpoint_residual;
        midpoint = 0.5 * (lower_bound + upper_bound);
        midpoint_residual = midpoint + benefit_cost + mu_rep *
          core_finite_mixture_moments(
            midpoint, u_sd, use_u_in_delta, component_mean,
            component_variance, component_weight
          )[3];
        if (lower_residual * midpoint_residual <= 0) {
          upper_bound = midpoint;
        } else {
          lower_bound = midpoint;
          lower_residual = midpoint_residual;
        }
      }
    } else {
      for (step in 1:120) {
        real residual = midpoint + benefit_cost + mu_rep *
          core_finite_mixture_moments(
            midpoint, u_sd, use_u_in_delta, component_mean,
            component_variance, component_weight
          )[3];
        midpoint = fmin(8, fmax(-8, midpoint - 0.4 * residual));
      }
    }
    return midpoint;
  }

  vector core_gq_calculate_roc(
      real benefit_cost,
      real benefit_cost_control,
      real total_error_sd,
      real u_sd,
      real dist_beta,
      real mu_rep,
      real mu_rep_control,
      real mu_rep_deriv,
      data int use_u_in_delta,
      data real alg_sol_rel_tol,
      data real alg_sol_f_tol,
      data real alg_sol_max_steps) {
    array[1] int x_i = {use_u_in_delta};
    real w = core_gq_find_fixedpoint_safe(
      benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta,
      alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps
    );
    real w_control = core_gq_find_fixedpoint_safe(
      benefit_cost_control, mu_rep_control, total_error_sd, u_sd,
      use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps
    );
    vector[2] delta_not_control = expected_delta_deriv(
      w, total_error_sd, u_sd, {0.0}, x_i
    );
    vector[2] delta_control = expected_delta_deriv(
      w_control, total_error_sd, u_sd, {0.0}, x_i
    );
    real roc = calculate_roc(
      w_control, total_error_sd, dist_beta, mu_rep, mu_rep_deriv,
      delta_control[1], delta_control[2]
    );
    real roc_no_visibility =
      -exp(normal_lpdf(w_control | 0, total_error_sd)) * dist_beta;
    real sm_value = calculate_sm(
      w, total_error_sd, dist_beta, mu_rep, mu_rep_deriv,
      delta_not_control[1], delta_not_control[2]
    );
    return [
      w, w_control, delta_control[1], delta_control[2], roc,
      roc_no_visibility, sm_value,
      delta_not_control[1], delta_not_control[2]
    ]';
  }

  vector core_gq_calculate_roc_student_t(
      real benefit_cost,
      real benefit_cost_control,
      real u_sd,
      real dist_beta,
      real mu_rep,
      real mu_rep_control,
      real mu_rep_deriv,
      data int use_u_in_delta,
      data real type_scale_sq,
      data vector type_precision,
      data vector type_weight) {
    real w = core_gq_find_fixedpoint_student_t_safe(
      benefit_cost, mu_rep, u_sd, use_u_in_delta, type_scale_sq,
      type_precision, type_weight
    );
    real w_control = core_gq_find_fixedpoint_student_t_safe(
      benefit_cost_control, mu_rep_control, u_sd, use_u_in_delta,
      type_scale_sq, type_precision, type_weight
    );
    vector[4] moments = core_student_t_moments(
      w, u_sd, use_u_in_delta, type_scale_sq,
      type_precision, type_weight
    );
    vector[4] control_moments = core_student_t_moments(
      w_control, u_sd, use_u_in_delta, type_scale_sq,
      type_precision, type_weight
    );
    real roc = -control_moments[2] *
      (dist_beta - mu_rep_deriv * control_moments[3]) /
      (1 + mu_rep * control_moments[4]);
    real roc_no_visibility = -control_moments[2] * dist_beta;
    real sm_value = calculate_sm(
      w, 1, dist_beta, mu_rep, mu_rep_deriv,
      moments[3], moments[4]
    );
    return [
      w, w_control, control_moments[3], control_moments[4], roc,
      roc_no_visibility, sm_value, moments[3], moments[4]
    ]';
  }

  vector core_gq_calculate_roc_finite_mixture(
      real benefit_cost, real benefit_cost_control, real u_sd,
      real dist_beta, real mu_rep, real mu_rep_control,
      real mu_rep_deriv, data int use_u_in_delta,
      vector component_mean, vector component_variance,
      vector component_weight) {
    real w = core_gq_find_fixedpoint_finite_mixture_safe(
      benefit_cost, mu_rep, u_sd, use_u_in_delta, component_mean,
      component_variance, component_weight
    );
    real w_control = core_gq_find_fixedpoint_finite_mixture_safe(
      benefit_cost_control, mu_rep_control, u_sd, use_u_in_delta,
      component_mean, component_variance, component_weight
    );
    vector[4] moments = core_finite_mixture_moments(
      w, u_sd, use_u_in_delta, component_mean,
      component_variance, component_weight
    );
    vector[4] control_moments = core_finite_mixture_moments(
      w_control, u_sd, use_u_in_delta, component_mean,
      component_variance, component_weight
    );
    real roc = -control_moments[2] *
      (dist_beta - mu_rep_deriv * control_moments[3]) /
      (1 + mu_rep * control_moments[4]);
    real roc_no_visibility = -control_moments[2] * dist_beta;
    real sm_value = calculate_sm(
      w, 1, dist_beta, mu_rep, mu_rep_deriv, moments[3], moments[4]
    );
    return [
      w, w_control, control_moments[3], control_moments[4], roc,
      roc_no_visibility, sm_value, moments[3], moments[4]
    ]';
  }

  real core_gq_calculate_noisy_sm(
      real benefit_cost,
      real lambda,
      real total_error_sd,
      real u_sd,
      real dist_beta,
      vector q_taker,
      vector q_nontaker,
      vector q_taker_dist,
      vector q_nontaker_dist,
      data int num_signals,
      data int use_u_in_delta,
      data real alg_sol_rel_tol,
      data real alg_sol_f_tol,
      data int alg_sol_max_steps) {
    array[1] int x_i = {use_u_in_delta};
    real cutoff = core_noisy_find_fixedpoint(
      benefit_cost, lambda, total_error_sd, u_sd, q_taker, q_nontaker,
      num_signals, use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol,
      alg_sol_max_steps
    );
    vector[2] delta = expected_delta_deriv(
      cutoff, total_error_sd, u_sd, {0.0}, x_i
    );
    vector[3] information = core_noisy_information_derivatives(
      cutoff, total_error_sd, q_taker, q_nontaker,
      q_taker_dist, q_nontaker_dist, num_signals
    );
    return (-dist_beta + lambda * information[3] * delta[1]) /
      (1 + lambda * (
        information[2] * delta[1] + information[1] * delta[2]
      ));
  }
}

data {
  #include struct_section_data.stan
  vector<lower=0>[num_clusters] core_cluster_weight;
  int<lower=0, upper=1> use_core_cluster_shock;
  real<lower=0> core_cluster_shock_sd_prior;
  int<lower=0, upper=2> core_lambda_structure;
  real<lower=0> core_lambda_log_ratio_sd_prior;
  int<lower=0, upper=1> core_profile_group_lambda;
  real core_profile_group_log_ratio;
  int<lower=0, upper=1> core_gq_override_lambda;
  vector<lower=0>[num_treatments] core_gq_lambda_override;
  int<lower=0, upper=2> core_type_distribution;
  real<lower=2> core_student_t_df;
  real<lower=0> core_type_scale_sq;
  int<lower=2> core_type_mixture_components;
  vector<lower=0>[core_type_mixture_components]
    core_type_mixture_precision;
  simplex[core_type_mixture_components] core_type_mixture_weight;
  array[num_wtp_obs] int<lower=1, upper=num_clusters> wtp_cluster_id;
  int<lower=0, upper=2> core_observation_model;
  int<lower=0, upper=2> core_recognition_structure;
  int<lower=0, upper=2> core_report_structure;
  int<lower=0, upper=1> core_report_arm_dist_hierarchical;
  real<lower=0> core_report_arm_dist_prior_scale;
  int<lower=0> core_num_peer_response_rows;
  array[core_num_peer_response_rows] int<lower=1, upper=num_clusters>
    core_peer_response_cluster_id;
  array[core_num_peer_response_rows] int<lower=0, upper=1>
    core_peer_true_takeup;
  array[core_num_peer_response_rows] int<lower=0> core_peer_total;
  array[core_num_peer_response_rows] int<lower=0> core_peer_recognized;
  array[core_num_peer_response_rows, 3] int<lower=0>
    core_peer_report_count;
}

transformed data {
  #include struct_section_transformed_data.stan
  array[num_clusters] int monitored_cluster_size = rep_array(0, num_clusters);
  array[num_clusters] int monitored_cluster_takeup = rep_array(0, num_clusters);
  array[num_treatments] int core_is_public_signal = rep_array(0, num_treatments);
  matrix[num_treatments, num_treatments - 1]
    core_signal_lambda_contrast_basis = rep_matrix(
      0,
      num_treatments,
      num_treatments - 1
    );
  for (contrast_index in 1:(num_treatments - 1)) {
    real denominator = sqrt(contrast_index * (contrast_index + 1.0));
    for (treatment_index in 1:contrast_index) {
      core_signal_lambda_contrast_basis[treatment_index, contrast_index] =
        1 / denominator;
    }
    core_signal_lambda_contrast_basis[contrast_index + 1, contrast_index] =
      -contrast_index / denominator;
  }
  core_is_public_signal[2] = 1;
  core_is_public_signal[BRACELET_TREATMENT_INDEX] = 1;
  if (num_treatments != 4 || CALENDAR_TREATMENT_INDEX != 3 ||
      BRACELET_TREATMENT_INDEX != 4) {
    reject("Core lambda grouping requires Control/Ink/Calendar/Bracelet order.");
  }
  if (core_observation_model > 0 &&
      (core_lambda_structure != 0 || core_type_distribution != 0)) {
    reject("Asymmetric observability requires common lambda and Gaussian types.");
  }
  if (core_observation_model == 0 &&
      (core_recognition_structure != 0 || core_report_structure != 0)) {
    reject("Observation restrictions require an asymmetric observation model.");
  }
  if (core_observation_model == 2 && core_recognition_structure == 2) {
    reject("The unconditional channel cannot condition recognition out.");
  }
  if (core_report_arm_dist_hierarchical &&
      (core_observation_model == 0 || core_report_structure != 0)) {
    reject("Hierarchical report-distance slopes require the full multinomial channel.");
  }
  for (included_index in 1:num_included_monitored_obs) {
    int obs_index = included_monitored_obs[included_index];
    int cluster_index = obs_cluster_id[obs_index];
    monitored_cluster_size[cluster_index] += 1;
    monitored_cluster_takeup[cluster_index] += takeup[obs_index];
  }
}

parameters {
  #include struct_section_parameters.stan
  vector[use_core_cluster_shock ? num_clusters : 0]
    core_cluster_shock_raw;
  vector<lower=0>[use_core_cluster_shock ? 1 : 0]
    core_cluster_shock_sd;
  vector[core_lambda_structure == 1 && !core_profile_group_lambda ? 1 : 0]
    core_lambda_group_log_ratio_raw;
  vector[core_lambda_structure == 2 ? num_treatments - 1 : 0]
    core_lambda_arm_log_ratio_raw;
  vector[core_observation_model > 0 && core_recognition_structure < 2 ? 2 : 0]
    core_recognition_intercept;
  vector[core_observation_model > 0 && core_recognition_structure == 0 ? 2 : 0]
    core_recognition_dist_slope;
  matrix[core_observation_model > 0 && core_recognition_structure == 0 ? 2 : 0, num_treatments - 1]
    core_recognition_arm_intercept_raw;
  matrix[core_observation_model > 0 && core_recognition_structure == 0 ? 2 : 0, num_treatments - 1]
    core_recognition_arm_dist_raw;
  matrix[core_observation_model > 0 && core_report_structure < 2 ? 2 : 0, 2]
    core_report_intercept;
  matrix[core_observation_model > 0 && core_report_structure < 2 ? 2 : 0, 2]
    core_report_dist_slope;
  matrix[core_observation_model > 0 && core_report_structure < 2 ? 2 : 0, 2 * (num_treatments - 1)]
    core_report_arm_intercept_raw;
  matrix[core_observation_model > 0 && core_report_structure == 0 ? 2 : 0, 2 * (num_treatments - 1)]
    core_report_arm_dist_raw;
  matrix<lower=0>[core_observation_model > 0 && core_report_structure == 0 && core_report_arm_dist_hierarchical ? 2 : 0, 2]
    core_report_arm_dist_sd;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0]
    core_definite_intercept;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0]
    core_definite_dist_slope;
  matrix[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0, num_treatments - 1]
    core_definite_arm_intercept_raw;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 1 : 0]
    core_definite_public_signal_dist_slope;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0]
    core_accuracy_intercept;
  matrix[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0, num_treatments - 1]
    core_accuracy_arm_intercept_raw;
  vector<lower=0, upper=1>[core_type_distribution == 2 ? 1 : 0]
    core_finite_mixture_weight;
  vector<lower=0, upper=1>[core_type_distribution == 2 ? 1 : 0]
    core_finite_mixture_between_share;
}

model {
}

generated quantities {
  // Rows: Combined, Close, Far. Columns follow the canonical Stan treatment
  // order and are reordered only by paper postprocessing.
  matrix[3, num_treatments] core_compact_takeup_level;
  // Rows: 500m, 1500m, 2500m.
  matrix[3, num_treatments] core_compact_sm;
  matrix[3, num_treatments] core_compact_sm_rescaled;
  // Secant multiplier from 500m to 2500m, normalized by the direct distance
  // cost change. This is reported beside the local derivative multiplier.
  vector[num_treatments] core_compact_sm_finite;
  // Rows are 500m, 1500m, 2500m; columns are treatment arms. Truth-specific
  // recognition and the equilibrium information factor make the fitted
  // response channel directly auditable without retaining latent arrays.
  matrix[3, num_treatments] core_compact_recognition_nontaker =
    rep_matrix(0, 3, num_treatments);
  matrix[3, num_treatments] core_compact_recognition_taker =
    rep_matrix(0, 3, num_treatments);
  matrix[3, num_treatments] core_compact_information_factor =
    rep_matrix(0, 3, num_treatments);
  matrix[3, num_treatments] core_compact_cutoff =
    rep_matrix(0, 3, num_treatments);
  matrix[3, num_treatments] core_compact_image_return =
    rep_matrix(0, 3, num_treatments);
  matrix[3, num_treatments] core_compact_sensitivity =
    rep_matrix(0, 3, num_treatments);
  matrix[3, num_treatments] core_compact_specificity =
    rep_matrix(0, 3, num_treatments);
  // Truth index 1=nontaker, 2=taker. Columns are Yes, No, Don't Know,
  // Unrecognized (the last is zero in the conditional-recognition mode).
  array[2] matrix[3 * num_treatments, 4] core_compact_response_matrix =
    rep_array(rep_matrix(0, 3 * num_treatments, 4), 2);
  vector[num_treatments] core_compact_signal_lambda =
    rep_vector(base_mu_rep, num_treatments);
  real core_compact_signal_vs_no_signal_log_lambda = 0;
  real core_compact_cluster_shock_sd = 0;
  real core_compact_log_lik_takeup = 0;
  real core_compact_log_lik_beliefs = 0;
  real core_compact_log_lik_wtp = 0;
  real core_compact_log_lik_distance = 0;

  {
    array[3] int compact_sm_roc_index = { 6, 16, 26 };
    array[num_discrete_dist] vector[num_dist_group_mix] group_dist_mean;
    vector[num_dist_group_treatments] beta =
      rep_vector(0, num_dist_group_treatments);
    vector[num_dist_group_treatments] structural_treatment_effect;
    vector[num_treatments] rep_intercept =
      beliefs_treatment_map_design_matrix * hyper_beta_1ord';
    vector[num_treatments] rep_dist_slope =
      beliefs_treatment_map_design_matrix * hyper_dist_beta_1ord';
    vector[num_clusters] cluster_shock = rep_vector(0, num_clusters);
    vector[num_treatments] combined_numerator =
      rep_vector(0, num_treatments);
    real total_weight = 0;
    real u_sd = raw_u_sd[1];
    real total_error_sd = sqrt(1 + square(u_sd));
    vector[2] finite_component_mean = rep_vector(0, 2);
    vector[2] finite_component_variance = rep_vector(1, 2);
    vector[2] finite_component_weight = rep_vector(0.5, 2);
    matrix[num_clusters, num_discrete_dist] all_cluster_dist;
    matrix[num_clusters, 2] observed_recognition = rep_matrix(0, num_clusters, 2);
    matrix[num_clusters, 3] observed_report_nontaker = rep_matrix(0, num_clusters, 3);
    matrix[num_clusters, 3] observed_report_taker = rep_matrix(0, num_clusters, 3);
    matrix[num_clusters, 4] observed_q_nontaker = rep_matrix(0, num_clusters, 4);
    matrix[num_clusters, 4] observed_q_taker = rep_matrix(0, num_clusters, 4);

    if (core_type_distribution == 2) {
      real mixture_weight = core_finite_mixture_weight[1];
      real between_share = core_finite_mixture_between_share[1];
      real separation = sqrt(
        between_share / (mixture_weight * (1 - mixture_weight))
      );
      finite_component_weight = [mixture_weight, 1 - mixture_weight]';
      finite_component_mean = [
        -(1 - mixture_weight) * separation,
        mixture_weight * separation
      ]';
      finite_component_variance = rep_vector(1 - between_share, 2);
    }

    if (core_lambda_structure == 1) {
      core_compact_signal_vs_no_signal_log_lambda =
        core_profile_group_lambda ? core_profile_group_log_ratio :
        core_lambda_log_ratio_sd_prior * core_lambda_group_log_ratio_raw[1];
      for (treatment_index in 1:num_treatments) {
        core_compact_signal_lambda[treatment_index] = base_mu_rep * exp(
          (core_is_public_signal[treatment_index] - 0.5) *
          core_compact_signal_vs_no_signal_log_lambda
        );
      }
    } else if (core_lambda_structure == 2) {
      core_compact_signal_lambda = base_mu_rep * exp(
        core_lambda_log_ratio_sd_prior / sqrt(2.0) *
        core_signal_lambda_contrast_basis *
        core_lambda_arm_log_ratio_raw
      );
      core_compact_signal_vs_no_signal_log_lambda = 0.5 * (
        log(core_compact_signal_lambda[2]) +
        log(core_compact_signal_lambda[BRACELET_TREATMENT_INDEX]) -
        log(core_compact_signal_lambda[1]) -
        log(core_compact_signal_lambda[CALENDAR_TREATMENT_INDEX])
      );
    }
    if (core_gq_override_lambda) {
      core_compact_signal_lambda = core_gq_lambda_override;
      if (min(core_compact_signal_lambda) > 0) {
        core_compact_signal_vs_no_signal_log_lambda = 0.5 * (
          log(core_compact_signal_lambda[2]) +
          log(core_compact_signal_lambda[BRACELET_TREATMENT_INDEX]) -
          log(core_compact_signal_lambda[1]) -
          log(core_compact_signal_lambda[CALENDAR_TREATMENT_INDEX])
        );
      } else {
        core_compact_signal_vs_no_signal_log_lambda = 0;
      }
    }

    if (use_core_cluster_shock) {
      core_compact_cluster_shock_sd = core_cluster_shock_sd[1];
      cluster_shock = core_cluster_shock_sd[1] * core_cluster_shock_raw;
    }

    for (dist_group_index in 1:num_discrete_dist) {
      group_dist_mean[dist_group_index] =
        hyper_dist_mean_mean +
        hyper_dist_mean_sd * group_dist_mean_raw[dist_group_index];
    }
    beta[1:2] = [beta_intercept, beta_ink_effect]';
    beta[CALENDAR_TREATMENT_INDEX] =
      beta_bracelet_effect + wtp_value_utility * hyper_wtp_mu;
    beta[BRACELET_TREATMENT_INDEX] = beta_bracelet_effect;
    structural_treatment_effect = treatment_map_design_matrix * beta;

    if (core_observation_model > 0) {
      for (cluster_index in 1:num_clusters) {
        int treatment_index = cluster_incentive_treatment_id[cluster_index];
        for (truth in 1:2) {
          vector[2] recognition = core_gq_recognition_channel(
            cluster_standard_dist[cluster_index], truth, treatment_index,
            core_recognition_structure, core_recognition_intercept,
            core_recognition_dist_slope, core_recognition_arm_intercept_raw,
            core_recognition_arm_dist_raw, core_signal_lambda_contrast_basis
          );
          vector[6] report = core_gq_report_channel(
            cluster_standard_dist[cluster_index], truth, treatment_index,
            core_report_structure, core_is_public_signal[treatment_index],
            core_report_intercept, core_report_dist_slope,
            core_report_arm_intercept_raw, core_report_arm_dist_raw,
            core_report_arm_dist_hierarchical, core_report_arm_dist_sd,
            core_definite_intercept, core_definite_dist_slope,
            core_definite_arm_intercept_raw,
            core_definite_public_signal_dist_slope, core_accuracy_intercept,
            core_accuracy_arm_intercept_raw,
            core_signal_lambda_contrast_basis
          );
          observed_recognition[cluster_index, truth] = recognition[1];
          if (truth == 1) {
            observed_report_nontaker[cluster_index] = report[1:3]';
            observed_q_nontaker[cluster_index] = core_noisy_signal_row(
              recognition[1], report[1:3], core_observation_model
            )';
          } else {
            observed_report_taker[cluster_index] = report[1:3]';
            observed_q_taker[cluster_index] = core_noisy_signal_row(
              recognition[1], report[1:3], core_observation_model
            )';
          }
        }
      }
    }

    if (fit_model_to_data) {
      vector[num_clusters] actual_mu =
        core_compact_signal_lambda[cluster_incentive_treatment_id] .* inv_logit(
          rep_intercept[cluster_incentive_treatment_id] +
          rep_dist_slope[cluster_incentive_treatment_id] .*
          cluster_standard_dist
        );
      vector[num_clusters] actual_benefit =
        structural_treatment_effect[cluster_assigned_dist_group_treatment] -
        dist_beta_v[1] * cluster_standard_dist + cluster_shock;
      for (cluster_index in included_clusters) {
        real cutoff;
        real takeup_probability;
        if (core_observation_model > 0) {
          cutoff = core_noisy_find_fixedpoint(
            actual_benefit[cluster_index], base_mu_rep, total_error_sd, u_sd,
            observed_q_taker[cluster_index]',
            observed_q_nontaker[cluster_index]',
            core_observation_model == 1 ? 3 : 4,
            use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol,
            alg_sol_max_steps
          );
          takeup_probability = Phi_approx(-cutoff / total_error_sd);
        } else if (core_type_distribution == 1) {
          cutoff = core_gq_find_fixedpoint_student_t_safe(
            actual_benefit[cluster_index], actual_mu[cluster_index], u_sd,
            use_u_in_delta, core_type_scale_sq,
            core_type_mixture_precision, core_type_mixture_weight
          );
          takeup_probability = 1 - core_student_t_moments(
            cutoff, u_sd, use_u_in_delta, core_type_scale_sq,
            core_type_mixture_precision, core_type_mixture_weight
          )[1];
        } else if (core_type_distribution == 2) {
          cutoff = core_gq_find_fixedpoint_finite_mixture_safe(
            actual_benefit[cluster_index], actual_mu[cluster_index], u_sd,
            use_u_in_delta, finite_component_mean,
            finite_component_variance, finite_component_weight
          );
          takeup_probability = 1 - core_finite_mixture_moments(
            cutoff, u_sd, use_u_in_delta, finite_component_mean,
            finite_component_variance, finite_component_weight
          )[1];
        } else {
          cutoff = core_gq_find_fixedpoint_safe(
            actual_benefit[cluster_index], actual_mu[cluster_index],
            total_error_sd, u_sd, use_u_in_delta, alg_sol_rel_tol,
            alg_sol_f_tol, alg_sol_max_steps
          );
          takeup_probability = Phi_approx(-cutoff / total_error_sd);
        }
        core_compact_log_lik_takeup += binomial_lpmf(
          monitored_cluster_takeup[cluster_index] |
          monitored_cluster_size[cluster_index],
          takeup_probability
        );
      }
    }

    if (fit_beliefs_model_to_data) {
      vector[num_beliefs_obs] belief_distance =
        cluster_standard_dist[beliefs_cluster_index];
      vector[num_beliefs_obs] predictor_1 =
        beliefs_treatment_design_matrix * hyper_beta_1ord' +
        (beliefs_treatment_design_matrix * hyper_dist_beta_1ord') .*
        belief_distance;
      vector[num_beliefs_obs] predictor_2 =
        beliefs_treatment_design_matrix * hyper_beta_2ord' +
        (beliefs_treatment_design_matrix * hyper_dist_beta_2ord') .*
        belief_distance;
      for (belief_index in 1:num_beliefs_obs) {
        if (core_observation_model == 0) {
          core_compact_log_lik_beliefs += binomial_logit_lpmf(
            num_knows_1ord[belief_index] | num_recognized[belief_index],
            predictor_1[belief_index]
          );
        }
        core_compact_log_lik_beliefs += binomial_logit_lpmf(
          num_knows_2ord[belief_index] | num_recognized[belief_index],
          predictor_2[belief_index]
        );
      }
      if (core_observation_model > 0) {
        for (row in 1:core_num_peer_response_rows) {
          int cluster_index = core_peer_response_cluster_id[row];
          int truth = core_peer_true_takeup[row] + 1;
          vector[3] report_probability = truth == 1 ?
            observed_report_nontaker[cluster_index]' :
            observed_report_taker[cluster_index]';
          core_compact_log_lik_beliefs += binomial_lpmf(
            core_peer_recognized[row] | core_peer_total[row],
            observed_recognition[cluster_index, truth]
          );
          core_compact_log_lik_beliefs += multinomial_lpmf(
            core_peer_report_count[row] | report_probability
          );
        }
      }
    }

    if (fit_dist_model_to_data) {
      for (cluster_index in 1:num_clusters) {
        int dist_group = cluster_treatment_map[
          cluster_assigned_dist_group_treatment[cluster_index], 2
        ];
        core_compact_log_lik_distance += lognormal_lpdf(
          cluster_standard_dist[cluster_index] |
          group_dist_mean[dist_group, 1], group_dist_sd[dist_group, 1]
        );
      }
    }

    if (fit_wtp_model_to_data) {
      int wtp_stratum_pos = 1;
      for (stratum_index in 1:num_strata) {
        int stratum_size = wtp_strata_sizes[stratum_index];
        int stratum_end = wtp_stratum_pos + stratum_size - 1;
        for (wtp_index in wtp_stratum_pos:stratum_end) {
          if (wtp_response[wtp_index] == -1) {
            if (gift_choice[wtp_index] == -1) {
              core_compact_log_lik_wtp += normal_lcdf(
                -scaled_wtp_offer[wtp_index] | hyper_wtp_mu, wtp_sigma
              );
            } else {
              core_compact_log_lik_wtp += normal_lccdf(
                scaled_wtp_offer[wtp_index] | hyper_wtp_mu, wtp_sigma
              );
            }
          } else {
            real endpoint = gift_choice[wtp_index] * scaled_wtp_offer[wtp_index];
            core_compact_log_lik_wtp += log_diff_exp(
              normal_lcdf(fmax(0, endpoint) | hyper_wtp_mu, wtp_sigma),
              normal_lcdf(fmin(0, endpoint) | hyper_wtp_mu, wtp_sigma)
            );
          }
        }
        wtp_stratum_pos = stratum_end + 1;
      }
    }

    // Match the paper's historical counterfactual-distance construction.
    for (cluster_index in 1:num_clusters) {
      int observed_dist_group = cluster_dist_treatment_id[cluster_index];
      real size_weight = sum(cluster_size[cluster_index]);
      total_weight += size_weight;
      for (dist_group_index in 1:num_discrete_dist) {
        if (dist_group_index == observed_dist_group) {
          all_cluster_dist[cluster_index, dist_group_index] =
            cluster_standard_dist[cluster_index];
        } else {
          all_cluster_dist[cluster_index, dist_group_index] = lognormal_rng(
            group_dist_mean[dist_group_index, 1],
            group_dist_sd[dist_group_index, 1]
          );
        }
      }
    }

    for (dist_group_index in 1:num_discrete_dist) {
      vector[num_clusters] counterfactual_dist =
        all_cluster_dist[, dist_group_index];
      for (treatment_index in 1:num_treatments) {
        int cell_index = 0;
        vector[num_clusters] mu_rep;
        vector[num_clusters] benefit_cost;
        vector[num_clusters] cutoff;
        vector[num_clusters] probability;
        matrix[num_clusters, 4] noisy_q_taker = rep_matrix(0, num_clusters, 4);
        matrix[num_clusters, 4] noisy_q_nontaker = rep_matrix(0, num_clusters, 4);
        real group_numerator = 0;

        for (candidate_index in 1:num_dist_group_treatments) {
          if (cluster_treatment_map[candidate_index, 1] == treatment_index &&
              cluster_treatment_map[candidate_index, 2] == dist_group_index) {
            cell_index = candidate_index;
          }
        }
        if (cell_index == 0) {
          reject("No treatment-by-distance cell found in compact core GQ.");
        }

        mu_rep = core_compact_signal_lambda[treatment_index] * inv_logit(
          rep_intercept[treatment_index] +
          rep_dist_slope[treatment_index] * counterfactual_dist
        );
        benefit_cost =
          rep_vector(structural_treatment_effect[cell_index], num_clusters) -
          dist_beta_v[1] * counterfactual_dist + cluster_shock;
        if (core_observation_model > 0) {
          for (cluster_index in 1:num_clusters) {
            for (truth in 1:2) {
              vector[2] recognition = core_gq_recognition_channel(
                counterfactual_dist[cluster_index], truth, treatment_index,
                core_recognition_structure, core_recognition_intercept,
                core_recognition_dist_slope,
                core_recognition_arm_intercept_raw,
                core_recognition_arm_dist_raw,
                core_signal_lambda_contrast_basis
              );
              vector[6] report = core_gq_report_channel(
                counterfactual_dist[cluster_index], truth, treatment_index,
                core_report_structure, core_is_public_signal[treatment_index],
                core_report_intercept, core_report_dist_slope,
                core_report_arm_intercept_raw, core_report_arm_dist_raw,
                core_report_arm_dist_hierarchical, core_report_arm_dist_sd,
                core_definite_intercept, core_definite_dist_slope,
                core_definite_arm_intercept_raw,
                core_definite_public_signal_dist_slope,
                core_accuracy_intercept, core_accuracy_arm_intercept_raw,
                core_signal_lambda_contrast_basis
              );
              if (truth == 1) {
                noisy_q_nontaker[cluster_index] = core_noisy_signal_row(
                  recognition[1], report[1:3],
                  core_observation_model
                )';
              } else {
                noisy_q_taker[cluster_index] = core_noisy_signal_row(
                  recognition[1], report[1:3],
                  core_observation_model
                )';
              }
            }
          }
        }
        for (cluster_index in 1:num_clusters) {
          if (core_observation_model > 0) {
            cutoff[cluster_index] = core_noisy_find_fixedpoint(
              benefit_cost[cluster_index], base_mu_rep, total_error_sd, u_sd,
              noisy_q_taker[cluster_index]', noisy_q_nontaker[cluster_index]',
              core_observation_model == 1 ? 3 : 4,
              use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol,
              alg_sol_max_steps
            );
          } else if (core_type_distribution == 1) {
            cutoff[cluster_index] = core_gq_find_fixedpoint_student_t_safe(
              benefit_cost[cluster_index], mu_rep[cluster_index], u_sd,
              use_u_in_delta, core_type_scale_sq,
              core_type_mixture_precision, core_type_mixture_weight
            );
          } else if (core_type_distribution == 2) {
            cutoff[cluster_index] =
              core_gq_find_fixedpoint_finite_mixture_safe(
                benefit_cost[cluster_index], mu_rep[cluster_index], u_sd,
                use_u_in_delta, finite_component_mean,
                finite_component_variance, finite_component_weight
              );
          } else {
            cutoff[cluster_index] = core_gq_find_fixedpoint_safe(
              benefit_cost[cluster_index],
              mu_rep[cluster_index],
              total_error_sd,
              u_sd,
              use_u_in_delta,
              alg_sol_rel_tol,
              alg_sol_f_tol,
              alg_sol_max_steps
            );
          }
        }
        if (core_type_distribution == 1) {
          for (cluster_index in 1:num_clusters) {
            probability[cluster_index] = 1 - core_student_t_moments(
              cutoff[cluster_index], u_sd, use_u_in_delta,
              core_type_scale_sq, core_type_mixture_precision,
              core_type_mixture_weight
            )[1];
          }
        } else if (core_type_distribution == 2) {
          for (cluster_index in 1:num_clusters) {
            probability[cluster_index] = 1 -
              core_finite_mixture_moments(
                cutoff[cluster_index], u_sd, use_u_in_delta,
                finite_component_mean, finite_component_variance,
                finite_component_weight
              )[1];
          }
        } else {
          probability = Phi_approx(-cutoff / total_error_sd);
        }

        for (cluster_index in 1:num_clusters) {
          real size_weight = sum(cluster_size[cluster_index]);
          group_numerator += size_weight * probability[cluster_index];
          if (cluster_dist_treatment_id[cluster_index] == dist_group_index) {
            combined_numerator[treatment_index] +=
              size_weight * probability[cluster_index];
          }
        }
        core_compact_takeup_level[dist_group_index + 1, treatment_index] =
          group_numerator / total_weight;
      }
    }
    core_compact_takeup_level[1] = combined_numerator' / total_weight;

    if (num_roc_distances < compact_sm_roc_index[3]) {
      reject("Compact core GQ requires ROC distances through 2500m.");
    }
    for (compact_dist_index in 1:3) {
      real roc_distance =
        roc_distances[compact_sm_roc_index[compact_dist_index]];
      vector[num_clusters] roc_cluster_dist =
        rep_vector(roc_distance, num_clusters);
      vector[num_clusters] control_mu;
      vector[num_clusters] control_benefit;
      int control_cell = 0;

      for (candidate_index in 1:num_dist_group_treatments) {
        if (cluster_treatment_map[candidate_index, 1] ==
            roc_compare_treatment_id_right) {
          control_cell = candidate_index;
        }
      }
      control_mu = core_compact_signal_lambda[
        roc_compare_treatment_id_right
      ] * inv_logit(
        rep_intercept[roc_compare_treatment_id_right] +
        rep_dist_slope[roc_compare_treatment_id_right] * roc_cluster_dist
      );
      control_benefit = rep_vector(
        structural_treatment_effect[control_cell], num_clusters
      ) - dist_beta_v[1] * roc_cluster_dist + cluster_shock;

      for (treatment_index in 1:num_treatments) {
        int treatment_cell = 0;
        vector[num_clusters] curr_mu;
        vector[num_clusters] curr_mu_deriv;
        vector[num_clusters] curr_mu_share;
        vector[num_clusters] curr_benefit;
        matrix[num_clusters, 9] roc_results;
        vector[4] noisy_q_taker = rep_vector(0, 4);
        vector[4] noisy_q_nontaker = rep_vector(0, 4);
        vector[4] noisy_q_taker_dist = rep_vector(0, 4);
        vector[4] noisy_q_nontaker_dist = rep_vector(0, 4);
        real noisy_information_sum = 0;
        real cutoff_sum = 0;
        real image_return_sum = 0;

        for (candidate_index in 1:num_dist_group_treatments) {
          if (cluster_treatment_map[candidate_index, 1] == treatment_index) {
            treatment_cell = candidate_index;
          }
        }
        curr_mu_share = inv_logit(
          rep_intercept[treatment_index] +
          rep_dist_slope[treatment_index] * roc_cluster_dist
        );
        curr_mu = core_compact_signal_lambda[treatment_index] * curr_mu_share;
        curr_mu_deriv = core_compact_signal_lambda[treatment_index] *
          rep_dist_slope[treatment_index] .*
          curr_mu_share .* (1 - curr_mu_share);
        curr_benefit = rep_vector(
          structural_treatment_effect[treatment_cell], num_clusters
        ) - dist_beta_v[1] * roc_cluster_dist + cluster_shock;

        if (core_observation_model > 0) {
          for (truth in 1:2) {
            vector[2] recognition = core_gq_recognition_channel(
              roc_distance, truth, treatment_index,
              core_recognition_structure, core_recognition_intercept,
              core_recognition_dist_slope,
              core_recognition_arm_intercept_raw,
              core_recognition_arm_dist_raw,
              core_signal_lambda_contrast_basis
            );
            vector[6] report = core_gq_report_channel(
              roc_distance, truth, treatment_index, core_report_structure,
              core_is_public_signal[treatment_index], core_report_intercept,
              core_report_dist_slope, core_report_arm_intercept_raw,
              core_report_arm_dist_raw, core_report_arm_dist_hierarchical,
              core_report_arm_dist_sd, core_definite_intercept,
              core_definite_dist_slope, core_definite_arm_intercept_raw,
              core_definite_public_signal_dist_slope,
              core_accuracy_intercept, core_accuracy_arm_intercept_raw,
              core_signal_lambda_contrast_basis
            );
            if (truth == 1) {
              core_compact_recognition_nontaker[
                compact_dist_index, treatment_index
              ] = recognition[1];
              if (core_observation_model == 1) {
                noisy_q_nontaker[1:3] = report[1:3];
                noisy_q_nontaker_dist[1:3] = report[4:6];
              } else {
                noisy_q_nontaker[1:3] =
                  recognition[1] * report[1:3];
                noisy_q_nontaker_dist[1:3] =
                  recognition[2] * report[1:3] +
                  recognition[1] * report[4:6];
                noisy_q_nontaker[4] = 1 - recognition[1];
                noisy_q_nontaker_dist[4] = -recognition[2];
              }
            } else {
              core_compact_recognition_taker[
                compact_dist_index, treatment_index
              ] = recognition[1];
              if (core_observation_model == 1) {
                noisy_q_taker[1:3] = report[1:3];
                noisy_q_taker_dist[1:3] = report[4:6];
              } else {
                noisy_q_taker[1:3] =
                  recognition[1] * report[1:3];
                noisy_q_taker_dist[1:3] =
                  recognition[2] * report[1:3] +
                  recognition[1] * report[4:6];
                noisy_q_taker[4] = 1 - recognition[1];
                noisy_q_taker_dist[4] = -recognition[2];
              }
            }
          }
          {
            int response_row =
              (compact_dist_index - 1) * num_treatments + treatment_index;
            core_compact_response_matrix[1][response_row] = noisy_q_nontaker';
            core_compact_response_matrix[2][response_row] = noisy_q_taker';
            core_compact_sensitivity[compact_dist_index, treatment_index] =
              noisy_q_taker[1];
            core_compact_specificity[compact_dist_index, treatment_index] =
              noisy_q_nontaker[2];
          }
        }

        for (cluster_index in 1:num_clusters) {
          if (core_observation_model > 0) {
            array[1] int delta_x_i = {use_u_in_delta};
            real cutoff = core_noisy_find_fixedpoint(
              curr_benefit[cluster_index], base_mu_rep,
              total_error_sd, u_sd, noisy_q_taker, noisy_q_nontaker,
              core_observation_model == 1 ? 3 : 4,
              use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol,
              alg_sol_max_steps
            );
            vector[2] delta = expected_delta_deriv(
              cutoff, total_error_sd, u_sd, {0.0}, delta_x_i
            );
            vector[3] information = core_noisy_information_derivatives(
              cutoff, total_error_sd, noisy_q_taker, noisy_q_nontaker,
              noisy_q_taker_dist, noisy_q_nontaker_dist,
              core_observation_model == 1 ? 3 : 4
            );
            roc_results[cluster_index] = rep_row_vector(0, 9);
            roc_results[cluster_index, 7] =
              (-dist_beta_v[1] + base_mu_rep * information[3] * delta[1]) /
              (1 + base_mu_rep * (
                information[2] * delta[1] + information[1] * delta[2]
              ));
            noisy_information_sum += information[1];
            cutoff_sum += cutoff;
            image_return_sum +=
              base_mu_rep * information[1] * delta[1];
          } else if (core_type_distribution == 1) {
            roc_results[cluster_index] = core_gq_calculate_roc_student_t(
              curr_benefit[cluster_index],
              control_benefit[cluster_index],
              u_sd,
              dist_beta_v[1],
              curr_mu[cluster_index],
              control_mu[cluster_index],
              curr_mu_deriv[cluster_index],
              use_u_in_delta,
              core_type_scale_sq,
              core_type_mixture_precision,
              core_type_mixture_weight
            )';
            cutoff_sum += roc_results[cluster_index, 1];
            image_return_sum +=
              curr_mu[cluster_index] * roc_results[cluster_index, 8];
          } else if (core_type_distribution == 2) {
            roc_results[cluster_index] =
              core_gq_calculate_roc_finite_mixture(
                curr_benefit[cluster_index], control_benefit[cluster_index],
                u_sd, dist_beta_v[1], curr_mu[cluster_index],
                control_mu[cluster_index], curr_mu_deriv[cluster_index],
                use_u_in_delta, finite_component_mean,
                finite_component_variance, finite_component_weight
              )';
            cutoff_sum += roc_results[cluster_index, 1];
            image_return_sum +=
              curr_mu[cluster_index] * roc_results[cluster_index, 8];
          } else {
            roc_results[cluster_index] = core_gq_calculate_roc(
              curr_benefit[cluster_index],
              control_benefit[cluster_index],
              total_error_sd,
              u_sd,
              dist_beta_v[1],
              curr_mu[cluster_index],
              control_mu[cluster_index],
              curr_mu_deriv[cluster_index],
              use_u_in_delta,
              alg_sol_rel_tol,
              alg_sol_f_tol,
              alg_sol_max_steps
            )';
            cutoff_sum += roc_results[cluster_index, 1];
            image_return_sum +=
              curr_mu[cluster_index] * roc_results[cluster_index, 8];
          }
        }
        core_compact_sm[compact_dist_index, treatment_index] =
          mean(roc_results[, 7]);
        core_compact_sm_rescaled[compact_dist_index, treatment_index] =
          core_compact_sm[compact_dist_index, treatment_index] /
          dist_beta_v[1];
        if (core_observation_model > 0) {
          core_compact_information_factor[compact_dist_index, treatment_index] =
            noisy_information_sum / num_clusters;
        }
        core_compact_cutoff[compact_dist_index, treatment_index] =
          cutoff_sum / num_clusters;
        core_compact_image_return[compact_dist_index, treatment_index] =
          image_return_sum / num_clusters;
      }
    }
    for (treatment_index in 1:num_treatments) {
      core_compact_sm_finite[treatment_index] =
        (core_compact_cutoff[3, treatment_index] -
         core_compact_cutoff[1, treatment_index]) /
        (dist_beta_v[1] *
         (roc_distances[compact_sm_roc_index[3]] -
          roc_distances[compact_sm_roc_index[1]]));
    }
  }
}
