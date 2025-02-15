row_vector[num_treatments] hyper_beta_1ord;
row_vector[beliefs_use_dist ? num_treatments : 0] hyper_dist_beta_1ord;
 
matrix[beliefs_use_stratum_level ? num_counties : 0, num_treatments] stratum_beta_1ord_raw;
matrix[beliefs_use_cluster_level ? num_clusters : 0, num_treatments] cluster_beta_1ord_raw;
matrix[beliefs_use_obs_level ? num_beliefs_obs : 0, num_treatments] obs_beta_1ord_raw;
 
row_vector<lower = 0>[beliefs_use_stratum_level ? num_treatments : 0] stratum_beta_1ord_sd;
row_vector<lower = 0>[beliefs_use_cluster_level ? num_treatments : 0] cluster_beta_1ord_sd;
row_vector<lower = 0>[beliefs_use_obs_level ? num_treatments : 0] obs_beta_1ord_sd;

matrix[beliefs_use_stratum_level && beliefs_use_dist ? num_counties : 0, num_treatments] stratum_dist_beta_1ord_raw;
matrix[beliefs_use_cluster_level && beliefs_use_dist ? num_clusters : 0, num_treatments] cluster_dist_beta_1ord_raw;
matrix[beliefs_use_obs_level ? num_beliefs_obs && beliefs_use_dist : 0, num_treatments] obs_dist_beta_1ord_raw;

row_vector<lower = 0>[beliefs_use_stratum_level && beliefs_use_dist ? num_treatments : 0] stratum_dist_beta_1ord_sd;
row_vector<lower = 0>[beliefs_use_cluster_level && beliefs_use_dist ? num_treatments : 0] cluster_dist_beta_1ord_sd;
row_vector<lower = 0>[beliefs_use_obs_level && beliefs_use_dist ? num_treatments : 0] obs_dist_beta_1ord_sd;
 
row_vector[num_treatments] hyper_beta_2ord;
row_vector[beliefs_use_dist ? num_treatments : 0] hyper_dist_beta_2ord;
 
matrix[beliefs_use_stratum_level? num_counties : 0, num_treatments] stratum_beta_2ord_raw;
matrix[beliefs_use_cluster_level ? num_clusters : 0, num_treatments] cluster_beta_2ord_raw;
matrix[beliefs_use_obs_level ? num_beliefs_obs : 0, num_treatments] obs_beta_2ord_raw;
 
row_vector<lower = 0>[beliefs_use_stratum_level ? num_treatments : 0] stratum_beta_2ord_sd;
row_vector<lower = 0>[beliefs_use_cluster_level ? num_treatments: 0] cluster_beta_2ord_sd;
row_vector<lower = 0>[beliefs_use_obs_level ? num_treatments : 0] obs_beta_2ord_sd;

matrix[beliefs_use_stratum_level && beliefs_use_dist ? num_counties : 0, num_treatments] stratum_dist_beta_2ord_raw;
matrix[beliefs_use_cluster_level && beliefs_use_dist ? num_clusters : 0, num_treatments] cluster_dist_beta_2ord_raw;
matrix[beliefs_use_obs_level && beliefs_use_dist ? num_beliefs_obs : 0, num_treatments] obs_dist_beta_2ord_raw;
 
row_vector<lower = 0>[beliefs_use_stratum_level && beliefs_use_dist ? num_treatments : 0] stratum_dist_beta_2ord_sd;
row_vector<lower = 0>[beliefs_use_cluster_level && beliefs_use_dist ? num_treatments: 0] cluster_dist_beta_2ord_sd;
row_vector<lower = 0>[beliefs_use_obs_level && beliefs_use_dist ? num_treatments : 0] obs_dist_beta_2ord_sd;
 
vector[beliefs_use_indiv_intercept ? num_beliefs_obs : 0] obs_beta_common_raw;
vector<lower = 0>[beliefs_use_indiv_intercept ? 1 : 0] obs_beta_common_sd;