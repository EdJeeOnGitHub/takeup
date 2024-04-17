
#include wtp_model_section.stan
#include beliefs_model_sec.stan
#include dist_model_sec.stan
#include struct_model_priors.stan


if (fit_model_to_data) {
// Take-up Likelihood 
if (use_binomial) {
    // Age groups not supported
    cluster_takeup_count[included_clusters, 1] ~ binomial(cluster_size[included_clusters, 1], structural_cluster_takeup_prob[included_clusters]);
} else {
    takeup[included_monitored_obs] ~ bernoulli(structural_cluster_takeup_prob[obs_cluster_id[included_monitored_obs]]);
}
}