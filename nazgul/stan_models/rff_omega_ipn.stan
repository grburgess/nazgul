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

  int<lower=1> k; // number of FFs

  int grainsize[N_detectors];
  int grainsize_meta;
  

}
transformed data {

  vector[max_N_time_bins] log_exposure[N_detectors] = log(exposure);
  
  real inv_sqrt_k = inv_sqrt(k);
  int detector_index[N_detectors];
  
  vector[3]  sc_diffs[N_detectors-1]; // precomputed sc position differences
  
  vector[N_detectors] maxs;
  vector[N_detectors] mins;

  // figure out the scales
  
  for (n in 1:N_detectors){

    maxs[n] = max(time[n]);
    mins[n] = min(time[n]);

    detector_index[n] = n;

  }

  for (n in 1:N_detectors){

    maxs[n] = max(time[n]);
    mins[n] = min(time[n]);

  }

  

  real max_range = max(maxs) - min(mins);


  // compute the differnces
  for (n in 1:N_detectors -1) {

    sc_diffs[n] = sc_pos[1] - sc_pos[n+1];


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

  for (n in 1:N_detectors-1) {

    
    dt[n] = time_delay(grb_xyz, sc_diffs[n]);

  }


}

model {

  vector[N_detectors] all_dt;
  vector[N_detectors] all_amplitude;

  all_dt[1] = 0.;
  all_amplitude[1] = 1.;

    
    
  
  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();

  /* log_scale[2] ~ normal(-1, 1.); */
  /* log_scale[1] ~ normal(0, 1); */

  raw_scale ~ normal(1,1);
  
  
  //log_bw ~ std_normal();

  //bw ~ cauchy(0, 2.5);

  range1_raw ~ lognormal(0, .2);
  //  range_delta ~ normal(0.5, 0.5);
  range2_raw ~ lognormal(log(1e-2), .2);

  //  range1_raw ~ exponential(1);
  //range2_raw ~ exponential(3);
  
  omega_var[1] ~ std_normal();
  omega_var[2] ~ std_normal();


  log_bkg ~ normal(log(500), log(100));  
  log_amplitude ~ std_normal();

  for (n in 2:N_detectors) {
    all_dt[n] = dt[n-1];
    all_amplitude[n] = amplitude[n-1];
  }



  
  /* target += reduce_sum(partial_log_like_bw_multi_scale_fast, counts[1,:N_time_bins[1]], grainsize[1], */
  /*                      time[1,:N_time_bins[1]], exposure[1,:N_time_bins[1]], */
  /*                      omega[1], omega[2], beta1, beta2, */
  /*                      0., bkg[1], scale[1], scale[2], 1, k); */


  /* target += reduce_sum(partial_log_like_bw_multi_scale_log, counts[1,:N_time_bins[1]], grainsize[1], */
  /*                      time[1,:N_time_bins[1]], log_exposure[1,:N_time_bins[1]], */
  /*                      omega[1], omega[2], beta1, beta2, */
  /*                      0., log_bkg[1], scale[1], scale[2], 0, k); */


  
  /* for (n in 2:N_detectors) { */

  /*   target += reduce_sum(partial_log_like_bw_multi_scale_fast, counts[n,:N_time_bins[n]], grainsize[n], */
  /*                        time[n,:N_time_bins[n]], exposure[n,:N_time_bins[n]], */
  /*                        omega[1], omega[2], beta1, beta2, */
  /*                        dt[n-1], bkg[n], scale[1], scale[2], amplitude[n-1], k); */

    /* target += reduce_sum(partial_log_like_bw_multi_scale_log, counts[n,:N_time_bins[n]], grainsize[n], */
    /*                      time[n,:N_time_bins[n]], log_exposure[n,:N_time_bins[n]], */
    /*                      omega[1], omega[2], beta1, beta2, */
    /*                      dt[n-1], log_bkg[n], scale[1], scale[2], log_amplitude[n-1], k); */

    
  /* } */


  target += reduce_sum(partial_total_like, detector_index, grainsize_meta,
		       counts, time, exposure,
		       omega[1], omega[2], beta1, beta2,
		       all_dt, bkg,
		       scale[1], scale[2],
		       all_amplitude, k, grainsize, N_time_bins);
  


}

generated quantities {

  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -( acos(grb_xyz[3]) - 0.5*pi());


}
