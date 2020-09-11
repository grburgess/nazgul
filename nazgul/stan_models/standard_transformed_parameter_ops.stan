
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

