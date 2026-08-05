vector core_finite_mixture_moments(
    real cutoff,
    real u_sd,
    int use_u_in_delta,
    vector component_mean,
    vector component_variance,
    vector component_weight) {
  int num_components = num_elements(component_weight);
  real below_probability = 0;
  real inference_numerator = 0;
  real index_density = 0;
  real inference_numerator_deriv = 0;
  real shock_variance = use_u_in_delta * square(u_sd);
  for (component in 1:num_components) {
    real type_mean = component_mean[component];
    real type_variance = component_variance[component];
    real index_variance = type_variance + shock_variance;
    real index_sd = sqrt(index_variance);
    real z = (cutoff - type_mean) / index_sd;
    real density_kernel = exp(std_normal_lpdf(z));
    real component_below = Phi_approx(z);
    below_probability += component_weight[component] * component_below;
    inference_numerator += component_weight[component] * (
      type_variance / index_sd * density_kernel -
      type_mean * component_below
    );
    index_density += component_weight[component] / index_sd * density_kernel;
    inference_numerator_deriv -= component_weight[component] *
      (type_variance * cutoff + type_mean * shock_variance) /
      (index_sd * index_variance) * density_kernel;
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

real core_finite_mixture_find_fixedpoint(
    real benefit_cost,
    real mu_rep,
    real u_sd,
    data int use_u_in_delta,
    vector component_mean,
    vector component_variance,
    vector component_weight) {
  real cutoff = fmin(8, fmax(-8, -benefit_cost));
  for (step_index in 1:10) {
    vector[4] moments = core_finite_mixture_moments(
      cutoff, u_sd, use_u_in_delta, component_mean,
      component_variance, component_weight
    );
    real residual = cutoff + benefit_cost + mu_rep * moments[3];
    real derivative = 1 + mu_rep * moments[4];
    real safe_derivative = abs(derivative) < 0.1 ?
      (derivative < 0 ? -0.1 : 0.1) : derivative;
    cutoff = fmin(8, fmax(-8, cutoff -
      fmin(1, fmax(-1, residual / safe_derivative))));
  }
  return cutoff;
}

vector core_finite_mixture_fixedpoint_rect(
    vector phi,
    vector theta,
    data array[] real x_r,
    data array[] int x_i) {
  int num_components = x_i[1];
  int use_u_in_delta = x_i[2];
  vector[num_components] component_mean;
  vector[num_components] component_variance;
  vector[num_components] component_weight;
  for (component in 1:num_components) {
    component_mean[component] = theta[3 + component];
    component_variance[component] = theta[3 + num_components + component];
    component_weight[component] = theta[3 + 2 * num_components + component];
  }
  return [core_finite_mixture_find_fixedpoint(
    theta[1], theta[2], theta[3], use_u_in_delta, component_mean,
    component_variance, component_weight
  )]';
}

vector core_finite_mixture_map_find_fixedpoint(
    vector benefit_cost,
    vector mu_rep,
    real u_sd,
    data int use_u_in_delta,
    vector component_mean,
    vector component_variance,
    vector component_weight) {
  int num_clusters = num_elements(benefit_cost);
  int num_components = num_elements(component_weight);
  vector[1] phi = [0.0]';
  array[num_clusters] vector[3 + 3 * num_components] theta;
  array[num_clusters, 2] int x_i;
  for (cluster in 1:num_clusters) {
    theta[cluster][1:3] = [benefit_cost[cluster], mu_rep[cluster], u_sd]';
    x_i[cluster, 1] = num_components;
    x_i[cluster, 2] = use_u_in_delta;
    for (component in 1:num_components) {
      theta[cluster][3 + component] = component_mean[component];
      theta[cluster][3 + num_components + component] =
        component_variance[component];
      theta[cluster][3 + 2 * num_components + component] =
        component_weight[component];
    }
  }
  return map_rect(
    core_finite_mixture_fixedpoint_rect, phi, theta,
    rep_array({0.0}, num_clusters), x_i
  );
}
