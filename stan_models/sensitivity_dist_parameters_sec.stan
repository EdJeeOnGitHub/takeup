// CmdStan CSV output is rounded, so a saved simplex can sum to
// 1 +/- a few 1e-7 and fail the strict simplex check during standalone
// generated quantities. For this read-only model, use a positive vector with
// the same names and dimensions. The original fitted values are otherwise
// unchanged.
array[num_discrete_dist] vector<lower = 0>[num_dist_group_mix] group_dist_mix;

array[num_discrete_dist] ordered[num_dist_group_mix] group_dist_mean_raw;
array[num_discrete_dist] vector<lower = 0>[num_dist_group_mix] group_dist_sd;

matrix[use_dist_county_effects ? num_discrete_dist : 0, num_counties]
  county_dist_effect_raw;
vector<lower = 0>[use_dist_county_effects ? num_discrete_dist : 0]
  county_dist_effect_sd;

matrix[use_dist_cluster_effects ? num_discrete_dist : 0, num_clusters]
  cluster_dist_effect_raw;
vector<lower = 0>[use_dist_county_effects ? num_discrete_dist : 0]
  cluster_dist_effect_sd;
