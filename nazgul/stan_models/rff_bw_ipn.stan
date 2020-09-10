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

  int grainsize;



}
transformed data {

  row_vector[k] omega1; // this weird MC integration thing.
  row_vector[k] omega2; // this weird MC integration thing.



  // for the non-delayed LC, let's go ahead and compute the fucking matrices

  for (i in 1:k) {

    omega1[i] = normal_rng(0, 1);
    omega2[i] = normal_rng(0, 1);
  }

}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  vector[N_detectors]  log_bkg;
  vector[N_detectors] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...

  real log_scale;
  
  real<lower=0> bw;

  unit_vector[3] grb_xyz; // GRB cartesian location


}

transformed parameters {

  vector[N_detectors] bkg = exp(log_bkg);
  vector[N_detectors] amplitude = exp(log_amplitude);

  real scale = exp(log_scale) * inv_sqrt(k);
 
  vector[N_detectors-1] dt;

  // compute all time delays relative to the first
  // detector

  for (n in 1:N_detectors-1) {

    
    dt[n] = time_delay(grb_xyz, sc_pos[1], sc_pos[n+1]);

  }


}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();

  log_scale ~ std_normal();

  log_bkg ~ normal(log(50), 1);
  bw ~ normal(5, 3);


  log_amplitude ~ std_normal();



  //  log_duration ~ normal(1,.2);
  //tstart ~ normal(1,5);

  target += reduce_sum(partial_log_like, counts[1], grainsize,
                       time[1], exposure[1],
                       omega1, omega2, beta1, beta2,bw,
                       0., bkg[1], scale, amplitude[1]);


  for (n in 2:N_detectors) {

    target += reduce_sum(partial_log_like, counts[n,:N_time_bins[n]], grainsize,
                         time[n,:N_time_bins[n]], exposure[n,:N_time_bins[n]],
                         omega1, omega2, beta1, beta2,bw,
                         dt[n-1], bkg[n], scale, amplitude[n]);

  }




}

generated quantities {

  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -( acos(grb_xyz[3]) - 0.5*pi());

  vector[k] omega[2];


  omega[1, :] = omega1';
  omega[2, :] = omega2';


}
