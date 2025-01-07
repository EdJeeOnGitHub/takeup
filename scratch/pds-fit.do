clear

import delimited "temp-data/analysis-cluster-covariate-data.csv"


list dewormed in 1/5

generate dewormed_num = (dewormed == "TRUE")
generate female_num = (female == "FALSE")
generate have_phone_num = (have_phone_lgl == "TRUE")
generate dist_pot_far = (distpotgroup == "far")


encode assignedtreatment, generate(assigned_treatment_fac)
recode assigned_treatment_fac (3 = 1 control) (4 = 2 ink) (2 = 3 calendar) (1 = 4 bracelet), gen(assigned_treatment)

encode county, generate(county_fac)





global cov_vars n_per_cluster      completed_primary         ///      
	 floor_tile_cement         ethnicity_luhya           ///
	 religion_christianity     treated_lgl               ///
	 adults_can_get_worms      correct_when_treat        ///
	 externality_omnibus       know_medicine_stops_worms ///
	 know_children_get_worms   sick_worms_only           ///
	 treated_past_year         praise_immuniz            ///
	 praise_dewor              stigma_immuniz            ///
	 stigma_dewor              female_num                ///
	 have_phone_num            agecensus                 


list $cov_vars in 1/5


help lincom


list county in 1/5

rename assigned_treatment a_treat
rename dist_pot_far dpf

reg dewormed_num i.a_treat##i.dpf $cov_vars i.county_fac, cluster(clusteridx)

ereturn list

matrix list e(b)

lincom _b[4.a_treat] - _b[4.a_treat#1.dpf]



reg dewormed_num i.assigned_treatment##c.standard_clusterdisttopot $cov_vars i.county_fac, cluster(clusteridx)

list clusterdisttopot




pdslasso dewormed_num  i.a_treat##i.dpf ($cov_vars i.county_fac), cluster(clusteridx) pnotpen(i.county_fac)

ereturn list

matrix list e(beta_plasso)


pdslasso dewormed_num  i.a_treat##c.standard_clusterdisttopot ($cov_vars i.county_fac), cluster(clusteridx) pnotpen(i.county_fac)

ssc install fre


fre distpotgroup

