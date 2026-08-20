// hierarchical linear regression w/ SS
functions {
    #include functions.stan
    real normal_ss(int N, real YY, row_vector YX, matrix XX, vector beta, real sigma) {
        real lp;
        lp = -0.5*N*log(sigma^2)  - YY * 0.5 * (1/sigma^2) + YX*beta*(1/sigma^2) - beta' * XX * beta * 0.5 * (1/sigma^2);
        return lp;
    }
}
data {
  int<lower=0> J; // number of groups
  int<lower=0> N; // number of observations total
  int<lower=1> K; // covariates/TEs
  
  vector[N] Y;
  matrix[N, K] X;
  array[N] int j_id;
  array[N] int signal_idx;
  array[4] int signal_x_idx;
  array[4] int no_signal_x_idx;


}
transformed data {
    array[J] real YY;
    array[J] row_vector[K] YX;
    array[J] matrix[K,K] XX;
    array[J] int N_j;
    array[J] int signal_idx_J;
    matrix[1, J] u = rep_matrix(1, 1, J);
    {
        matrix[N, 1] Y_mat = to_matrix(Y);
        for (site_j in 1:J) {
            N_j[site_j] = num_test(j_id, { site_j }, 1);
            array[N_j[site_j]] int j_ids = which(j_id, { site_j }, 1);
            YY[site_j] = (Y_mat[j_ids, ]' * Y_mat[j_ids, ])[1,1];
            YX[site_j] = to_row_vector(Y_mat[j_ids, ]' * X[j_ids, ]);
            XX[site_j] = X[j_ids, ]' * X[j_ids,];
            signal_idx_J[site_j] = signal_idx[j_ids][1];
        }
    }
}

parameters {
  vector<lower=0>[K] sigma_beta;                       // prediction error scale
  array[J] real<lower=0> sigma_y;
  matrix[J, K] z;
  vector[K] gamma;
  vector[4] mu_signal;
  vector[4] mu_no_signal;

}
transformed parameters {
    matrix[J, K] beta;
    for (j in 1:J) {
        beta[j, ] = to_row_vector(gamma + z[j , ]' .* sigma_beta); 
    }
}
model {
    
  // priors
  mu_no_signal ~ normal(0, 1);
  mu_signal ~ normal(0, 1);
  gamma[signal_x_idx] ~ normal(mu_signal, 1);
  gamma[no_signal_x_idx] ~ normal(mu_no_signal, 1);
  to_vector(z) ~ std_normal();
  sigma_beta ~ std_normal();
  sigma_y ~ std_normal();

  for (j in 1:J){
    target += normal_ss(N_j[j], YY[j], YX[j], XX[j], to_vector(beta[j]), sigma_y[J]);
  }  
}
generated quantities {
}
