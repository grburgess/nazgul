functions {
#include functions.stan
}

data {


  int N_detectors; // number of detectors used
  int N_time_bins[N_detectors]; // the number of time bins in each detector
  int max_N_time_bins; // the max number of time bins

  int counts[N_detectors, max_N_time_bins]; // matrix of counts per detector

  vector[max_N_time_bins] time[N_detectors]; // matrix of time bin midpoints per detector
  vector[max_N_time_bins] exposure[N_detectors]; // matrix of exposures per detector


  vector[3] sc_pos[N_detectors]; // the 3D vector for each space craft
  vector[3] sc_pointing_norm[N_detectors]; // the 3D vector for each space craft

  int<lower=1> k; // number of FFs

  int grainsize[N_detectors];



}
transformed data {

  vector[max_N_time_bins] log_exposure[N_detectors] = log(exposure);
  
  real inv_sqrt_k = inv_sqrt(k);

  vector[N_detectors] maxs;
  vector[N_detectors] mins;

  for (n in 1:N_detectors){

    maxs[n] = max(time[n]);
    mins[n] = min(time[n]);

  }


  real max_range = max(maxs) - min(mins);

  vector[N_detectors] horizon_angle;
  vector[3] sc_pos_norm[N_detectors];

  for (n in 1:N_detectors) {

    horizon_angle[n] = calculate_horizon_angle(sc_pos[n]);
    sc_pos_norm[n] = norm_vector(sc_pos[n]);
    
  }

  
}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  row_vector[k] omega_var[2]; // this weird MC integration thing.
  



  vector[N_detectors]  log_bkg;
  vector[N_detectors-1] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...
  
  positive_ordered [2] raw_scale;
  //real<lower=0> log_scale;
  //positive_ordered[2] bw;
  //ordered[2] log_bw;

  real<lower=0, upper=1> range1_raw;
  real<lower=0, upper=range1_raw> range2_raw;
  

  unit_vector[3] grb_xyz; // GRB cartesian location


}

transformed parameters {

  vector[N_detectors] bkg = exp(log_bkg);
  vector[N_detectors-1] amplitude = exp(log_amplitude);
  vector[N_detectors] amplitude_mod;
  vector[N_detectors] earth_occulted;;
  
  
   vector[2] scale = raw_scale * inv_sqrt(k);
  //real scale = log_scale * inv_sqrt(k);

  vector[2] range;
  vector[2] bw;
  

  row_vector[k] omega[2]; // this weird MC integration thing.
  
  vector[N_detectors-1] dt;

  range[1] = range1_raw * max_range;
  range[2] = range2_raw * max_range;

  bw = inv(range);

  
  // non-center
  omega[1] = omega_var[1] * bw[1];
  omega[2] = omega_var[2] * bw[2];
  
  // compute all time delays relative to the first
  // detector

  for (n in 1:N_detectors) {

    earth_occulted[n]= earth_occulation(horizon_angle[n], sc_pos_norm[n], grb_xyz);
    
  }
  
  for (n in 1:N_detectors-1) {

    
    dt[n] = time_delay(grb_xyz, sc_pos[1], sc_pos[n+1]);

  }

  amplitude_mod[1] = earth_occulted[1];

  for (n in 1:N_detectors- 1){

    amplitude_mod[n+1] = amplitude[n] * earth_occulted[n];
  }

}

model {

  
  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();

  /* log_scale[2] ~ normal(-1, 1.); */
  /* log_scale[1] ~ normal(0, 1); */

  raw_scale ~ normal(1,1);
  
  
  //log_bw ~ std_normal();

  //bw ~ cauchy(0, 2.5);

  range1_raw ~ lognormal(log(1e-1), .5);
  //  range_delta ~ normal(0.5, 0.5);
  range2_raw ~ lognormal(log(1e-1), .5);

  
  omega_var[1] ~ std_normal();
  omega_var[2] ~ std_normal();


  log_bkg ~ normal(log(500), log(100));  
  log_amplitude ~ std_normal();





  


  
  target += reduce_sum(partial_log_like_bw_multi_scale_fast, counts[1,:N_time_bins[1]], grainsize[1],
                       time[1,:N_time_bins[1]], exposure[1,:N_time_bins[1]],
                       omega[1], omega[2], beta1, beta2,
                       0., bkg[1], scale[1], scale[2], amplitude_mod[1], k);


  /* target += reduce_sum(partial_log_like_bw_multi_scale_log, counts[1,:N_time_bins[1]], grainsize[1], */
  /*                      time[1,:N_time_bins[1]], log_exposure[1,:N_time_bins[1]], */
  /*                      omega[1], omega[2], beta1, beta2, */
  /*                      0., log_bkg[1], scale[1], scale[2], 0, k); */


  
  for (n in 2:N_detectors) {

    target += reduce_sum(partial_log_like_bw_multi_scale_fast, counts[n,:N_time_bins[n]], grainsize[n],
                         time[n,:N_time_bins[n]], exposure[n,:N_time_bins[n]],
                         omega[1], omega[2], beta1, beta2,
                         dt[n-1], bkg[n], scale[1], scale[2], amplitude_mod[n], k);

    /* target += reduce_sum(partial_log_like_bw_multi_scale_log, counts[n,:N_time_bins[n]], grainsize[n], */
    /*                      time[n,:N_time_bins[n]], log_exposure[n,:N_time_bins[n]], */
    /*                      omega[1], omega[2], beta1, beta2, */
    /*                      dt[n-1], log_bkg[n], scale[1], scale[2], log_amplitude[n-1], k); */

    
  }




}

generated quantities {

  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -( acos(grb_xyz[3]) - 0.5*pi());


}
