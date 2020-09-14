
  vector[max_N_time_bins] log_exposure[N_detectors] = log(exposure);
  
  real inv_sqrt_k = inv_sqrt(k);
  int detector_index[N_detectors];
  
  vector[3]  sc_diffs[N_detectors-1]; // precomputed sc position differences
  
  vector[N_detectors] maxs;
  vector[N_detectors] mins;


