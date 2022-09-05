functions {
  // Differences between inverse Gamma tail probabilities and target probabilities
  vector tail_delta(vector y, vector theta, array[] real x_r, array[] int x_i) {
    vector[2] deltas;
    deltas[1] = inv_gamma_cdf(theta[1] | exp(y[1]), exp(y[2])) - x_r[1];
    deltas[2] = 1 - inv_gamma_cdf(theta[2] | exp(y[1]), exp(y[2])) - x_r[2];
    return deltas;
  }
}

transformed data {
  // Target quantiles
  real l = 0.125; // Lower quantile
  real u = 2;   // Upper quantile
  vector[2] theta = [l, u]';
  
  // Initial guess at inverse Gamma parameters 
  // using asymmetric Gaussian approximation
  real dl = 0.2;
  real du = 5;
  real alpha_guess = square(2 * (dl * u + du * l) / (u - l)) + 2;
  real beta_guess = (alpha_guess - 1) * 0.5 * (dl * u + du * l) / (dl + du);
  vector[2] y_guess = [log(alpha_guess), log(beta_guess)]';

  // Find inverse Gamma density parameters that ensure 
  // 1% probabilty below l and 1% probabilty above u
  vector[2] y;
  array[2] real x_r = { 0.01, 0.01 };
  array[0] int x_i;

  y = algebra_solver(tail_delta, y_guess, theta, x_r, x_i);

  print("alpha = ", exp(y[1]));
  print("beta = ", exp(y[2]));
}

generated quantities {
  real alpha = exp(y[1]);
  real beta = exp(y[2]);
}