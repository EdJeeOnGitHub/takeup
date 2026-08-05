real core_noisy_information_factor(
    real cutoff,
    real total_error_sd,
    vector q_taker,
    vector q_nontaker,
    data int num_signals) {
  real takeup_prob = Phi_approx(-cutoff / total_error_sd);
  real factor_sum = 0;
  for (signal in 1:num_signals) {
    real difference = q_taker[signal] - q_nontaker[signal];
    real signal_prob = takeup_prob * q_taker[signal] +
      (1 - takeup_prob) * q_nontaker[signal];
    factor_sum += square(difference) / signal_prob;
  }
  return takeup_prob * (1 - takeup_prob) * factor_sum;
}

// Returns A, dA/dw, and dA/ddistance for the binary latent action observed
// through an arbitrary discrete response channel. q_*_dist are the analytic
// distance derivatives of the response probabilities.
vector core_noisy_information_derivatives(
    real cutoff,
    real total_error_sd,
    vector q_taker,
    vector q_nontaker,
    vector q_taker_dist,
    vector q_nontaker_dist,
    data int num_signals) {
  real z = -cutoff / total_error_sd;
  real takeup_prob = Phi_approx(z);
  real takeup_prob_w = -takeup_prob * (1 - takeup_prob) *
    (1.5976 + 3 * 0.07056 * square(z)) / total_error_sd;
  real h = 0;
  real h_prob = 0;
  real h_dist = 0;
  for (signal in 1:num_signals) {
    real difference = q_taker[signal] - q_nontaker[signal];
    real difference_dist =
      q_taker_dist[signal] - q_nontaker_dist[signal];
    real signal_prob = takeup_prob * q_taker[signal] +
      (1 - takeup_prob) * q_nontaker[signal];
    real signal_prob_dist = takeup_prob * q_taker_dist[signal] +
      (1 - takeup_prob) * q_nontaker_dist[signal];
    h += square(difference) / signal_prob;
    h_prob -= difference ^ 3 / square(signal_prob);
    h_dist += 2 * difference * difference_dist / signal_prob -
      square(difference) * signal_prob_dist / square(signal_prob);
  }
  return [
    takeup_prob * (1 - takeup_prob) * h,
    takeup_prob_w * (
      (1 - 2 * takeup_prob) * h +
      takeup_prob * (1 - takeup_prob) * h_prob
    ),
    takeup_prob * (1 - takeup_prob) * h_dist
  ]';
}

vector core_noisy_fixedpoint_residual(
    vector cutoff,
    vector theta,
    data array[] real x_r,
    data array[] int x_i) {
  int use_u_in_delta = x_i[1];
  int num_signals = x_i[2];
  vector[4] q_taker = theta[5:8];
  vector[4] q_nontaker = theta[9:12];
  real delta = expected_delta(
    cutoff[1], theta[3], theta[4], x_r, {use_u_in_delta}
  );
  real information_factor = core_noisy_information_factor(
    cutoff[1], theta[3], q_taker, q_nontaker, num_signals
  );
  return [cutoff[1] + theta[1] +
    theta[2] * information_factor * delta]';
}

real core_noisy_find_fixedpoint(
    real benefit_cost,
    real lambda,
    real total_error_sd,
    real u_sd,
    vector q_taker,
    vector q_nontaker,
    data int num_signals,
    data int use_u_in_delta,
    data real alg_sol_rel_tol,
    data real alg_sol_f_tol,
    data int alg_sol_max_steps) {
  vector[12] theta;
  theta[1:4] = [benefit_cost, lambda, total_error_sd, u_sd]';
  theta[5:8] = q_taker;
  theta[9:12] = q_nontaker;
  return algebra_solver(
    core_noisy_fixedpoint_residual,
    [-benefit_cost]',
    theta,
    {0.0},
    {use_u_in_delta, num_signals},
    alg_sol_rel_tol,
    alg_sol_f_tol,
    alg_sol_max_steps
  )[1];
}

vector core_noisy_find_fixedpoint_rect(
    vector phi,
    vector theta,
    data array[] real x_r,
    data array[] int x_i) {
  return [core_noisy_find_fixedpoint(
    theta[1], theta[2], theta[3], theta[4], theta[5:8], theta[9:12],
    x_i[2], x_i[1], x_r[1], x_r[2], x_i[3]
  )]';
}

vector core_noisy_map_find_fixedpoint(
    vector benefit_cost,
    real lambda,
    vector total_error_sd,
    vector u_sd,
    matrix q_taker,
    matrix q_nontaker,
    data int num_signals,
    data int use_u_in_delta,
    data real alg_sol_rel_tol,
    data real alg_sol_f_tol,
    data int alg_sol_max_steps) {
  int num_clusters = num_elements(benefit_cost);
  array[num_clusters] vector[12] theta;
  array[num_clusters, 3] int x_i = rep_array(
    {use_u_in_delta, num_signals, alg_sol_max_steps}, num_clusters
  );
  for (cluster in 1:num_clusters) {
    theta[cluster][1:4] = [
      benefit_cost[cluster], lambda, total_error_sd[cluster], u_sd[cluster]
    ]';
    theta[cluster][5:8] = q_taker[cluster]';
    theta[cluster][9:12] = q_nontaker[cluster]';
  }
  return map_rect(
    core_noisy_find_fixedpoint_rect,
    [total_error_sd[1]]',
    theta,
    rep_array({alg_sol_rel_tol, alg_sol_f_tol}, num_clusters),
    x_i
  );
}

vector core_softmax_with_reference(real yes_logit, real no_logit) {
  return softmax([yes_logit, no_logit, 0]');
}

vector core_noisy_signal_row(
    real recognition_prob,
    vector conditional_report,
    data int observation_model) {
  vector[4] result = rep_vector(0, 4);
  if (observation_model == 1) {
    result[1:3] = conditional_report;
  } else {
    result[1:3] = recognition_prob * conditional_report;
    result[4] = 1 - recognition_prob;
  }
  return result;
}

vector core_two_stage_report_row(
    real definite_prob,
    real accuracy_prob,
    data int truth) {
  vector[3] result;
  if (truth == 2) {
    result = [
      definite_prob * accuracy_prob,
      definite_prob * (1 - accuracy_prob),
      1 - definite_prob
    ]';
  } else {
    result = [
      definite_prob * (1 - accuracy_prob),
      definite_prob * accuracy_prob,
      1 - definite_prob
    ]';
  }
  return result;
}

vector core_two_stage_report_row_derivative(
    real definite_prob,
    real definite_slope,
    real accuracy_prob,
    data int truth) {
  real definite_deriv = definite_prob * (1 - definite_prob) * definite_slope;
  vector[3] result;
  if (truth == 2) {
    result = [
      definite_deriv * accuracy_prob,
      definite_deriv * (1 - accuracy_prob),
      -definite_deriv
    ]';
  } else {
    result = [
      definite_deriv * (1 - accuracy_prob),
      definite_deriv * accuracy_prob,
      -definite_deriv
    ]';
  }
  return result;
}
