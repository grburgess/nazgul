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


  int<lower=1> k; // number of FFs

  int grainsize[N_detectors];


}
transformed data {

  vector[max_N_time_bins] log_exposure[N_detectors] = log(exposure);
  
  /* vector[2] bw; */
  /* bw[1] = 1.5; */
  /* bw[2] = 0.3; */

  real inv_sqrt_k = inv_sqrt(k);

  vector[N_detectors] maxs;
  vector[N_detectors] mins;

  for (n in 1:N_detectors){

    maxs[n] = max(time[n]);
    mins[n] = min(time[n]);

  }


  real max_range = max(maxs) - min(mins);

  
  
}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  row_vector[k] omega_var[2]; // this weird MC integration thing.

  real  log_bkg;
  
  vector[2] log_scale;
  //vector<lower=0>[2] scale_raw;
  
  //positive_ordered[2] bw;
  //ordered[2] log_bw;

  real<lower=0, upper=1> range1_raw;
  real<lower=0, upper=1> range2_raw;
  //real<lower=0,upper=1> range_delta;
  
}

transformed parameters {

  //  real bkg = exp(log_bkg);

  vector[2] scale = exp(log_scale) * inv_sqrt_k;
  //vector[2] scale = scale_raw * inv_sqrt(k);

  vector[2] range;

  vector[2] bw;
  
  // vector[2] bw = exp(log_bw);
  
  row_vector[k] omega[2]; // this weird MC integration thing.

  range[1] = range1_raw * max_range;
  //range[2] = range[1] * range_delta;
  range[2] = range2_raw * max_range;

  bw = inv(range);

  
  // non-center
  omega[1] = omega_var[1] * bw[1];
  omega[2] = omega_var[2] * bw[2];
  
  
}

model {

  // priors

  
  beta1 ~ std_normal();
  beta2 ~ std_normal();

  log_scale[2] ~ normal(-2, .5);
  log_scale[1] ~ normal(-1, .5);
  


  /* log_bw[1] ~ normal(-1, .5); */
  /* log_bw[2] ~ normal(0, .5); */

  range1_raw ~ normal(0, 1);
  //  range_delta ~ normal(0.5, 0.5);
  range2_raw ~ normal(0, 1);
  

  //bw ~ std_normal();
  
  omega_var[1] ~ std_normal();
  omega_var[2] ~ std_normal();

  log_bkg ~ normal(log(500), log(100));
  
  /* target += reduce_sum(partial_log_like_bw_multi_scale, counts[1,:N_time_bins[1]], grainsize[1], */
  /*                      time[1,:N_time_bins[1]], exposure[1,:N_time_bins[1]], */
  /*                      omega[1], omega[2], beta1, beta2, */
  /*                      0., bkg, scale[1], scale[2], 1., k); */



  target += reduce_sum(partial_log_like_bw_multi_scale_log, counts[1,:N_time_bins[1]], grainsize[1],
                       time[1,:N_time_bins[1]], log_exposure[1,:N_time_bins[1]],
                       omega[1], omega[2], beta1, beta2,
                       0., log_bkg, scale[1], scale[2], 0, k);

  
}

generated quantities {

  real bkg = exp(log_bkg);

  
}
