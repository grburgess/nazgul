all_amplitude[1] = ang_sep[1] * earth_occulted[1];
for (n in 2:N_detectors) {
  all_dt[n] = dt[n-1];
  all_amplitude[n] = amplitude[n-1] * ang_sep[n-1] * earth_occulted[n-1];
 }
