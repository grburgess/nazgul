
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

