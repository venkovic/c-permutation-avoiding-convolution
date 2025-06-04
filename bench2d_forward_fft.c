/*
 * Copyright (c) 2025 Nicolas Venkovic
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

  int t2_1 = 12;           // Default: 12
  int timing_runs = 10;    // Default: 10

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
      t2_1 = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
      timing_runs = atoi(argv[++i]);
    } else {
      fprintf(stderr, "Unknown or malformed option: %s\n", argv[i]);
      return 1;
    }
  }

  int t2_2 = t2_1;

  int n1 = 1 << t2_1; // n1 = 2^t2_1
  int n2 = 1 << t2_2; // n2 = 2^t2_2
  int n  = n1 * n2;

  // Allocate memory
  double* re  = malloc(n * sizeof(double));
  double* im  = malloc(n * sizeof(double));
  fftw_complex* fftw_in  = fftw_malloc(n * sizeof(fftw_complex));
  fftw_complex* fftw_out = fftw_malloc(n * sizeof(fftw_complex));
  double* fftw_in_flat  = (double*)fftw_in;
  double* fftw_out_flat = (double*)fftw_out;

  int* rho = precompute_index_reversal_permutation_r2_2d(n1, n2);
  double** twiddle_re1, **twiddle_im1;
  double** twiddle_re2, **twiddle_im2;
  precompute_twiddles_r2(n1, &twiddle_re1, &twiddle_im1);
  precompute_twiddles_r2(n2, &twiddle_re2, &twiddle_im2);
  fftw_plan p = fftw_plan_dft_2d(n1, n2, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);

  if (do_verify) {
    initialize_data_split(re, im, fftw_in_flat, n);
    
    // FFTW
    fftw_execute(p);

    // Radix-2
    fft_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho);

    double max_err = get_max_error(re, im, fftw_out_flat, n);
    printf("Max error FFT radix-2 vs FFTW: %.3e\n", max_err);

    free(rho);
    free_twiddles(twiddle_re1, twiddle_im1, t2_1);
    free_twiddles(twiddle_re2, twiddle_im2, t2_2);

    // Radix-4
    if (supports_radix_4(t2_1) && supports_radix_4(t2_2)) {
      int t4_1 = t2_1 >> 1;
      int t4_2 = t2_2 >> 1;
      int* rho = precompute_index_reversal_permutation_r4_2d(n1, n2);

      double** twiddle_re1, **twiddle_im1;
      double** twiddle_re2, **twiddle_im2;
      precompute_twiddles_r4(n1, &twiddle_re1, &twiddle_im1);
      precompute_twiddles_r4(n2, &twiddle_re2, &twiddle_im2);
          
      initialize_data_split(re, im, fftw_in_flat, n);
        
      fft_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho);
      fftw_execute(p);

      max_err = get_max_error(re, im, fftw_out_flat, n);
      printf("Max error FFT radix-4 vs FFTW: %.3e\n", max_err);          

      free(rho);
      free_twiddles(twiddle_re1, twiddle_im1, t4_1);
      free_twiddles(twiddle_re2, twiddle_im2, t4_2);
    }
  }

  if (do_bench) {
    struct timespec start, end;

    // Radix-2
    
    // Warm-up
    initialize_data_split(re, im, fftw_in_flat, n);
    fft_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho);
    apply_permutation_split(re, im, n, rho);
    ufft_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2);
    
    // Timing runs
    double total_fft = 0., total_perm = 0., total_ufft = 0.;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fft_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_fft += time_diff(start, end);

      clock_gettime(CLOCK_MONOTONIC, &start);
      apply_permutation_split(re, im, n, rho);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_perm += time_diff(start, end);

      clock_gettime(CLOCK_MONOTONIC, &start);
      ufft_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_ufft += time_diff(start, end);
    }

    printf("Average radix-2 FFT time:            %.6f seconds\n", total_fft / timing_runs);
    printf("Average radix-2 permutation time:    %.6f seconds\n", total_perm / timing_runs);
    printf("Average radix-2 unordered FFT time:  %.6f seconds\n", total_ufft / timing_runs);
    
    // Clean-up
    free(rho);
    free_twiddles(twiddle_re1, twiddle_im1, t2_1);  
    free_twiddles(twiddle_re2, twiddle_im2, t2_2);  

    // Radix-4
    if (supports_radix_4(t2_1) && supports_radix_4(t2_2)) {
      int t4_1 = t2_1 >> 1;
      int t4_2 = t2_2 >> 1;
      int* rho = precompute_index_reversal_permutation_r4_2d(n1, n2);
      double** twiddle_re1, **twiddle_im1;
      double** twiddle_re2, **twiddle_im2;
      precompute_twiddles_r4(n, &twiddle_re1, &twiddle_im1);
      precompute_twiddles_r4(n, &twiddle_re2, &twiddle_im2);

      // Warm-up
      initialize_data_split(re, im, fftw_in_flat, n);
      fft_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho);
      apply_permutation_split(re, im, n, rho);
      ufft_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2);

      // Timing runs
      total_fft = 0., total_perm = 0., total_ufft = 0.;
      for (int t = 0; t < timing_runs; ++t) {
        initialize_data_split(re, im, fftw_in_flat, n);

        clock_gettime(CLOCK_MONOTONIC, &start);
        fft_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_fft += time_diff(start, end);

        clock_gettime(CLOCK_MONOTONIC, &start);
        apply_permutation_split(re, im, n, rho);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_perm += time_diff(start, end);

        clock_gettime(CLOCK_MONOTONIC, &start);
        ufft_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_ufft += time_diff(start, end);
      }

      printf("Average radix-4 FFT time:            %.6f seconds\n", total_fft / timing_runs);
      printf("Average radix-4 permutation time:    %.6f seconds\n", total_perm / timing_runs);
      printf("Average radix-4 unordered FFT time:  %.6f seconds\n", total_ufft / timing_runs);
        
      // Clean-up
      free(rho);
      free_twiddles(twiddle_re1, twiddle_im1, t4_1); 
      free_twiddles(twiddle_re2, twiddle_im2, t4_2); 
    }

    // FFTW

    // Warm-up
    initialize_data_split(re, im, fftw_in_flat, n);
    fftw_execute(p);

    // Timing runs
    double total_fftw = 0.;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fftw_execute(p);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_fftw += time_diff(start, end);
    }

    printf("Average FFTW time:                   %.6f seconds\n", total_fftw / timing_runs);
  }

  // Clean-up
  fftw_destroy_plan(p);
  fftw_free(fftw_in);
  fftw_free(fftw_out);
  free(re); free(im); 

  return 0;
}