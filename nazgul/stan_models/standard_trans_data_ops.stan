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
