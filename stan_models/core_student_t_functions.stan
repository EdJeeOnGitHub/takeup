vector core_student_t_moments(
    real cutoff,
    real u_sd,
    int use_u_in_delta,
    real type_scale_sq,
    vector type_precision,
    vector type_weight) {
  int num_components = num_elements(type_weight);
  real below_probability = 0;
  real inference_numerator = 0;
  real index_density = 0;
  real inference_numerator_deriv = 0;

  // A variance-one t_nu type is a normal variance mixture:
  // v | q ~ N(0, ((nu - 2) / nu) / q), q ~ Gamma(nu/2, nu/2).
  // The supplied generalized Gauss-Laguerre rule integrates over q.
  for (component in 1:num_components) {
    real type_var = type_scale_sq / type_precision[component];
    real index_var = type_var + use_u_in_delta * square(u_sd);
    real index_sd = sqrt(index_var);
    real standardized_cutoff = cutoff / index_sd;
    real density_kernel = exp(std_normal_lpdf(standardized_cutoff));
    below_probability += type_weight[component] *
      Phi_approx(standardized_cutoff);
    inference_numerator += type_weight[component] *
      type_var / index_sd * density_kernel;
    index_density += type_weight[component] / index_sd * density_kernel;
    inference_numerator_deriv -= type_weight[component] *
      type_var * cutoff / (index_sd * index_var) * density_kernel;
  }
  below_probability = fmin(1 - 1e-12, fmax(1e-12, below_probability));
  {
    real denominator = below_probability * (1 - below_probability);
    real delta = inference_numerator / denominator;
    real delta_deriv = inference_numerator_deriv / denominator -
      inference_numerator * index_density * (1 - 2 * below_probability) /
      square(denominator);
    return [below_probability, index_density, delta, delta_deriv]';
  }
}

vector core_student_t_fixedpoint_system(
    vector cutoff,
    vector theta,
    data array[] real x_r,
    data array[] int x_i) {
  int num_components = x_i[1];
  int use_u_in_delta = x_i[2];
  vector[num_components] precision;
  vector[num_components] weight;
  for (component in 1:num_components) {
    precision[component] = theta[4 + component];
    weight[component] = theta[4 + num_components + component];
  }
  {
    vector[4] moments = core_student_t_moments(
      cutoff[1], theta[3], use_u_in_delta, theta[4], precision, weight
    );
    return [cutoff[1] + theta[1] + theta[2] * moments[3]]';
  }
}

real core_student_t_find_fixedpoint(
    real benefit_cost,
    real mu_rep,
    real u_sd,
    data int use_u_in_delta,
    data real type_scale_sq,
    data vector type_precision,
    data vector type_weight,
    data real alg_sol_rel_tol,
    data real alg_sol_f_tol,
    data real alg_sol_max_steps) {
  real cutoff = fmin(8, fmax(-8, -benefit_cost));
  for (step_index in 1:8) {
    vector[4] moments = core_student_t_moments(
      cutoff, u_sd, use_u_in_delta, type_scale_sq,
      type_precision, type_weight
    );
    real residual = cutoff + benefit_cost + mu_rep * moments[3];
    real residual_deriv = 1 + mu_rep * moments[4];
    real safe_deriv = abs(residual_deriv) < 0.1 ?
      (residual_deriv < 0 ? -0.1 : 0.1) : residual_deriv;
    real newton_step = fmin(1, fmax(-1, residual / safe_deriv));
    cutoff = fmin(8, fmax(-8, cutoff - newton_step));
  }
  return cutoff;
}

vector core_student_t_fixedpoint_rect(
    vector phi,
    vector theta,
    data array[] real x_r,
    data array[] int x_i) {
  int num_components = x_i[1];
  int use_u_in_delta = x_i[2];
  vector[num_components] precision;
  vector[num_components] weight;
  real cutoff = fmin(8, fmax(-8, -theta[1]));
  for (component in 1:num_components) {
    precision[component] = theta[4 + component];
    weight[component] = theta[4 + num_components + component];
  }
  for (step_index in 1:8) {
    vector[4] moments = core_student_t_moments(
      cutoff, theta[3], use_u_in_delta, theta[4], precision, weight
    );
    real residual = cutoff + theta[1] + theta[2] * moments[3];
    real residual_deriv = 1 + theta[2] * moments[4];
    real safe_deriv = abs(residual_deriv) < 0.1 ?
      (residual_deriv < 0 ? -0.1 : 0.1) : residual_deriv;
    real newton_step = fmin(1, fmax(-1, residual / safe_deriv));
    cutoff = fmin(8, fmax(-8, cutoff - newton_step));
  }
  return [cutoff]';
}

vector core_student_t_map_find_fixedpoint(
    vector benefit_cost,
    vector mu_rep,
    real u_sd,
    data int use_u_in_delta,
    data real type_scale_sq,
    data vector type_precision,
    data vector type_weight,
    data real alg_sol_rel_tol,
    data real alg_sol_f_tol,
    data real alg_sol_max_steps) {
  int num_clusters = num_elements(benefit_cost);
  int num_components = num_elements(type_weight);
  vector[1] phi = [0.0]';
  array[num_clusters] vector[4 + 2 * num_components] theta;
  array[num_clusters, 2] int x_i;
  for (cluster in 1:num_clusters) {
    theta[cluster][1:4] = [
      benefit_cost[cluster], mu_rep[cluster], u_sd, type_scale_sq
    ]';
    x_i[cluster, 1] = num_components;
    x_i[cluster, 2] = use_u_in_delta;
    for (component in 1:num_components) {
      theta[cluster][4 + component] = type_precision[component];
      theta[cluster][4 + num_components + component] = type_weight[component];
    }
  }
  return map_rect(
    core_student_t_fixedpoint_rect, phi, theta,
    rep_array({0.0}, num_clusters), x_i
  );
}
