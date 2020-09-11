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
