/*
 * Copyright (c) 2026 Nicolas Venkovic
 * 
 * This file is part of c-permutation-avoiding-convolution.
 * 
 * This file is licensed under the MIT License.
 * For the full license text, see the LICENSE file in the root directory of this project.
 */

#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <fftw3.h>

#include "fft.h"
           
int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s [verify|bench] [-t power] [-r timing_runs]\n", argv[0]);
    return 1;
  }

  int do_verify = strcmp(argv[1], "verify") == 0;
  int do_bench = strcmp(argv[1], "bench") == 0;

  if (!do_verify && !do_bench) {
    fprintf(stderr, "Unknown mode: %s\n", argv[1]);
    return 1;
  }

  int t2 = 24;           // Default: 24
  int timing_runs = 10;  // Default: 10

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
      t2 = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
      timing_runs = atoi(argv[++i]);
    } else {
      fprintf(stderr, "Unknown or malformed option: %s\n", argv[i]);
      return 1;
    }
  }

  int n = 1 << t2; // n = 2^t2

  // Allocate memory
  double* re  = malloc(n * sizeof(double));
  double* im  = malloc(n * sizeof(double));
  fftw_complex* fftw_in  = fftw_malloc(n * sizeof(fftw_complex));
  fftw_complex* fftw_out = fftw_malloc(n * sizeof(fftw_complex));
  double* fftw_in_flat  = (double*)fftw_in;
  double* fftw_out_flat = (double*)fftw_out;

  int* rho = precompute_index_reversal_permutation_r2(n);
  double** twiddle_re, **twiddle_im;
  precompute_twiddles_r2(n, &twiddle_re, &twiddle_im);
  fftw_plan p = fftw_plan_dft_1d(n, fftw_in, fftw_out, FFTW_BACKWARD, FFTW_ESTIMATE);

  if (do_verify) {
    initialize_data_split(re, im, fftw_in_flat, n);
    
    // FFTW
    fftw_execute(p);
    for (int i = 0; i < n; i++) {
      fftw_out_flat[2 * i]     /= n;
      fftw_out_flat[2 * i + 1] /= n;
    }

    // Radix-2
    ifft(re, im, n, 2, twiddle_re, twiddle_im, rho);

    double max_err = get_max_error(re, im, fftw_out_flat, n);
    printf("Max error FFT radix-2 vs FFTW: %.3e\n", max_err);

    free(rho);
    free_twiddles(twiddle_re, twiddle_im, t2);

    // Radix-4
    if (supports_radix_4(t2)) {
      int t4 = t2 >> 1;
      int* rho = precompute_index_reversal_permutation_r4(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r4(n, &twiddle_re, &twiddle_im);
          
      initialize_data_split(re, im, fftw_in_flat, n);
        
      ifft(re, im, n, 4, twiddle_re, twiddle_im, rho);
      fftw_execute(p);
      for (int i = 0; i < n; i++) {
        fftw_out_flat[2 * i]     /= n;
        fftw_out_flat[2 * i + 1] /= n;
      }

      max_err = get_max_error(re, im, fftw_out_flat, n);
      printf("Max error FFT radix-4 vs FFTW: %.3e\n", max_err);          

      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t4);
    }

    // Radix-8
    if (supports_radix_8(t2)) {
      int t8 = t2 >> 2;
      int* rho = precompute_index_reversal_permutation_r8(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r8(n, &twiddle_re, &twiddle_im);

      initialize_data_split(re, im, fftw_in_flat, n);
        
      ifft(re, im, n, 8, twiddle_re, twiddle_im, rho);
      fftw_execute(p);
      for (int i = 0; i < n; i++) {
        fftw_out_flat[2 * i]     /= n;
        fftw_out_flat[2 * i + 1] /= n;
      }

      max_err = get_max_error(re, im, fftw_out_flat, n);
      printf("Max error FFT radix-8 vs FFTW: %.3e\n", max_err);

      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t8);      
    }
  }

  if (do_bench) {
    struct timespec start, end;

    // Radix-2
    
    // Warm-up
    initialize_data_split(re, im, fftw_in_flat, n);
    ifft(re, im, n, 2, twiddle_re, twiddle_im, rho);
    apply_permutation_split(re, im, n, rho);
    uifft(re, im, n, 2, twiddle_re, twiddle_im);
    
    // Timing runs
    double total_fft = 0., total_perm = 0., total_ufft = 0.;
    double total_fft_sq = 0., total_perm_sq = 0., total_ufft_sq = 0.;
    double dt;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      ifft(re, im, n, 2, twiddle_re, twiddle_im, rho);
      clock_gettime(CLOCK_MONOTONIC, &end);
      dt = time_diff(start, end);
      total_fft += dt;
      total_fft_sq += dt * dt;

      clock_gettime(CLOCK_MONOTONIC, &start);
      apply_permutation_split(re, im, n, rho);
      clock_gettime(CLOCK_MONOTONIC, &end);
      dt = time_diff(start, end);
      total_perm += dt;
      total_perm_sq += dt * dt;

      clock_gettime(CLOCK_MONOTONIC, &start);
      uifft(re, im, n, 2, twiddle_re, twiddle_im);
      clock_gettime(CLOCK_MONOTONIC, &end);
      dt = time_diff(start, end);
      total_ufft += dt;
      total_ufft_sq += dt * dt;
    }

    double mu = total_fft / timing_runs;
    double s = sqrt(total_fft_sq / timing_runs - mu * mu);
    printf("Average radix-2 FFT time:             %.6f +/- %.6f seconds\n", mu, s);

    mu = total_perm / timing_runs;
    s = sqrt(total_perm_sq / timing_runs - mu * mu);
    printf("Average radix-2 permutation time:     %.6f +/- %.6f seconds\n", mu, s);

    mu = total_ufft / timing_runs;
    s = sqrt(total_ufft_sq / timing_runs - mu * mu);
    printf("Average radix-2 unordered FFT time:   %.6f +/- %.6f seconds\n", mu, s);
    
    // Clean-up
    free(rho);
    free_twiddles(twiddle_re, twiddle_im, t2);  

    // Radix-4
    if (supports_radix_4(t2)) {
      int t4 = t2 >> 1;
      int* rho = precompute_index_reversal_permutation_r4(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r4(n, &twiddle_re, &twiddle_im);
          
      // Warm-up
      initialize_data_split(re, im, fftw_in_flat, n);
      ifft(re, im, n, 4, twiddle_re, twiddle_im, rho);
      apply_permutation_split(re, im, n, rho);
      uifft(re, im, n, 4, twiddle_re, twiddle_im);

      // Timing runs
      total_fft = 0., total_perm = 0., total_ufft = 0.;
      total_fft_sq = 0., total_perm_sq = 0., total_ufft_sq = 0.;
      for (int t = 0; t < timing_runs; ++t) {
        initialize_data_split(re, im, fftw_in_flat, n);

        clock_gettime(CLOCK_MONOTONIC, &start);
        ifft(re, im, n, 4, twiddle_re, twiddle_im, rho);
        clock_gettime(CLOCK_MONOTONIC, &end);
        dt = time_diff(start, end);
        total_fft += dt;
        total_fft_sq += dt * dt;

        clock_gettime(CLOCK_MONOTONIC, &start);
        apply_permutation_split(re, im, n, rho);
        clock_gettime(CLOCK_MONOTONIC, &end);  
        dt = time_diff(start, end);
        total_perm += dt;
        total_perm_sq += dt * dt;

        clock_gettime(CLOCK_MONOTONIC, &start);
        uifft(re, im, n, 4, twiddle_re, twiddle_im);
        clock_gettime(CLOCK_MONOTONIC, &end);
        dt = time_diff(start, end);
        total_ufft += dt;
        total_ufft_sq += dt * dt;
      }

      mu = total_fft / timing_runs;
      s = sqrt(total_fft_sq / timing_runs - mu * mu);
      printf("Average radix-4 FFT time:             %.6f +/- %.6f seconds\n", mu, s);

      mu = total_perm / timing_runs;
      s = sqrt(total_perm_sq / timing_runs - mu * mu);
      printf("Average radix-4 permutation time:     %.6f +/- %.6f seconds\n", mu, s);

      mu = total_ufft / timing_runs;
      s = sqrt(total_ufft_sq / timing_runs - mu * mu);
      printf("Average radix-4 unordered FFT time:   %.6f +/- %.6f seconds\n", mu, s);
        
      // Clean-up
      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t4); 
    }

    // Radix-8
    if (supports_radix_8(t2)) {
      int t8 = t2 >> 2;
      int* rho = precompute_index_reversal_permutation_r8(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r8(n, &twiddle_re, &twiddle_im);
          
      // Warm-up
      initialize_data_split(re, im, fftw_in_flat, n);
      ifft(re, im, n, 8, twiddle_re, twiddle_im, rho);
      apply_permutation_split(re, im, n, rho);
      uifft(re, im, n, 8, twiddle_re, twiddle_im);

      // Timing runs
      total_fft = 0., total_perm = 0., total_ufft = 0.;
      total_fft_sq = 0., total_perm_sq = 0., total_ufft_sq = 0.;
      for (int t = 0; t < timing_runs; ++t) {
        initialize_data_split(re, im, fftw_in_flat, n);

        clock_gettime(CLOCK_MONOTONIC, &start);
        ifft(re, im, n, 8, twiddle_re, twiddle_im, rho);
        clock_gettime(CLOCK_MONOTONIC, &end);
        dt = time_diff(start, end);
        total_fft += dt;
        total_fft_sq += dt * dt;

        clock_gettime(CLOCK_MONOTONIC, &start);
        apply_permutation_split(re, im, n, rho);
        clock_gettime(CLOCK_MONOTONIC, &end);
        dt = time_diff(start, end);
        total_perm += dt;
        total_perm_sq += dt * dt;

        clock_gettime(CLOCK_MONOTONIC, &start);
        uifft(re, im, n, 8, twiddle_re, twiddle_im);
        clock_gettime(CLOCK_MONOTONIC, &end);
        dt = time_diff(start, end);
        total_ufft += dt;
        total_ufft_sq += dt * dt;
      }

      mu = total_fft / timing_runs;
      s = sqrt(total_fft_sq / timing_runs - mu * mu);
      printf("Average radix-8 FFT time:             %.6f +/- %.6f seconds\n", mu, s);

      mu = total_perm / timing_runs;
      s = sqrt(total_perm_sq / timing_runs - mu * mu);
      printf("Average radix-8 permutation time:     %.6f +/- %.6f seconds\n", mu, s);

      mu = total_ufft / timing_runs;
      s = sqrt(total_ufft_sq / timing_runs - mu * mu);
      printf("Average radix-8 unordered FFT time:   %.6f +/- %.6f seconds\n", mu, s);

      // Clean-up
      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t8); 
    }

    // FFTW

    // Warm-up
    initialize_data_split(re, im, fftw_in_flat, n);
    fftw_execute(p);
    for (int i = 0; i < n; i++) {
      fftw_out_flat[2 * i]     /= n;
      fftw_out_flat[2 * i + 1] /= n;
    }

    // Timing runs
    double total_fftw = 0., total_fftw_sq = 0.;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fftw_execute(p);
      for (int i = 0; i < n; i++) {
        fftw_out_flat[2 * i]     /= n;
        fftw_out_flat[2 * i + 1] /= n;
      }
      clock_gettime(CLOCK_MONOTONIC, &end);
      dt = time_diff(start, end);
      total_fftw += dt;
      total_fftw_sq += dt * dt;
    }

    mu = total_fftw / timing_runs;
    s = sqrt(total_fftw_sq / timing_runs - mu * mu);
    printf("Average FFTW time:                    %.6f +/- %.6f seconds\n", mu, s);
  }

  // Clean-up
  fftw_destroy_plan(p);
  fftw_free(fftw_in);
  fftw_free(fftw_out);
  free(re); free(im); 

  return 0;
}
