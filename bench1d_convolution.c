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

void fftw_based_convolution(fftw_plan p_forward, fftw_plan p_backward,
                            double* fftw_in_flat, double* fftw_out_flat,
                            double* h_re, double* h_im,
                            double* re, double* im,
                            int n) {
  
  fftw_execute(p_forward);

  for (int i = 0; i < n; i++) {
    re[i] = fftw_out_flat[2 * i];
    im[i] = fftw_out_flat[2 * i + 1];

    fftw_in_flat[2 * i] = h_re[i];
    fftw_in_flat[2 * i + 1] = h_im[i];
  }

  fftw_execute(p_forward);

  for (int i = 0; i < n; i++) {
    double Fx_re = re[i];
    double Fx_im = im[i];

    double Fh_re = fftw_out_flat[2 * i];
    double Fh_im = fftw_out_flat[2 * i + 1];

    fftw_out_flat[2 * i]     = (Fx_re * Fh_re - Fx_im * Fh_im);
    fftw_out_flat[2 * i + 1] = (Fx_re * Fh_im + Fx_im * Fh_re);
  }

  fftw_execute(p_backward);

  for (int i = 0; i < n; i++) {
    fftw_in_flat[2 * i] /= n;
    fftw_in_flat[2 * i + 1] /= n;
  }
}
           
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
  double* tmp_re  = malloc(n * sizeof(double));
  double* tmp_im  = malloc(n * sizeof(double));
  double* h_re  = malloc(n * sizeof(double));
  double* h_im  = malloc(n * sizeof(double));
  double* h_re_copy  = malloc(n * sizeof(double));
  double* h_im_copy  = malloc(n * sizeof(double));
  fftw_complex* fftw_in  = fftw_malloc(n * sizeof(fftw_complex));
  fftw_complex* fftw_out = fftw_malloc(n * sizeof(fftw_complex));
  double* fftw_in_flat  = (double*)fftw_in;
  double* fftw_out_flat = (double*)fftw_out;

  int* rho = precompute_index_reversal_permutation_r2(n);
  double** twiddle_re, **twiddle_im;
  precompute_twiddles_r2(n, &twiddle_re, &twiddle_im);

  fftw_plan p_forward  = fftw_plan_dft_1d(n, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan p_backward = fftw_plan_dft_1d(n, fftw_out, fftw_in, FFTW_BACKWARD, FFTW_ESTIMATE);

  initialize_filter(h_re, h_im, n);
  memcpy(h_re_copy, h_re, n * sizeof(double));
  memcpy(h_im_copy, h_im, n * sizeof(double));

  if (do_verify) {
    initialize_data_split(re, im, fftw_in_flat, n);
    
    // FFTW
    fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);

    // FFT-based radix-2
    fft_based_convolution(re, im, n, 2, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);

    double max_err = get_max_error(re, im, fftw_in_flat, n);
    printf("Max error FFT-based radix-2 vs FFTW: %.3e\n", max_err);

    initialize_data_split(re, im, fftw_in_flat, n);
    
    // FFTW
    fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);

    memcpy(h_re_copy, h_re, n * sizeof(double));
    memcpy(h_im_copy, h_im, n * sizeof(double));

    // Permutation-avoiding radix-2
    permutation_avoiding_convolution(re, im, n, 2, twiddle_re, twiddle_im, h_re_copy, h_im_copy);

    max_err = get_max_error(re, im, fftw_in_flat, n);
    printf("Max error Permutation-avoiding radix-2 vs FFTW: %.3e\n", max_err);

    free(rho);
    free_twiddles(twiddle_re, twiddle_im, t2);

    // Radix-4
    if (supports_radix_4(t2)) {
      int t4 = t2 >> 1;
      int* rho = precompute_index_reversal_permutation_r4(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r4(n, &twiddle_re, &twiddle_im);
      
      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));

      initialize_data_split(re, im, fftw_in_flat, n);
      
      // FFTW
      fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);

      // FFT-based radix-4
      fft_based_convolution(re, im, n, 4, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);

      max_err = get_max_error(re, im, fftw_in_flat, n);
      printf("Max error FFT-based radix-4 vs FFTW: %.3e\n", max_err);          

      initialize_data_split(re, im, fftw_in_flat, n);
    
      // FFTW
      fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);

      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));

      // Permutation-avoiding radix-4
      permutation_avoiding_convolution(re, im, n, 4, twiddle_re, twiddle_im, h_re_copy, h_im_copy);

      max_err = get_max_error(re, im, fftw_in_flat, n);
      printf("Max error Permutation-avoiding radix-4 vs FFTW: %.3e\n", max_err);

      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t4);
    }

    // Radix-8
    if (supports_radix_8(t2)) {
      int t8 = t2 >> 2;
      int* rho = precompute_index_reversal_permutation_r8(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r8(n, &twiddle_re, &twiddle_im);
      
      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));

      initialize_data_split(re, im, fftw_in_flat, n);
      
      // FFTW
      fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);

      // FFT-based radix-8
      fft_based_convolution(re, im, n, 8, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);

      max_err = get_max_error(re, im, fftw_in_flat, n);
      printf("Max error FFT-based radix-8 vs FFTW: %.3e\n", max_err);          

      initialize_data_split(re, im, fftw_in_flat, n);
    
      // FFTW
      fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);

      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));

      // Permutation-avoiding radix-8
      permutation_avoiding_convolution(re, im, n, 8, twiddle_re, twiddle_im, h_re_copy, h_im_copy);

      max_err = get_max_error(re, im, fftw_in_flat, n);
      printf("Max error Permutation-avoiding radix-8 vs FFTW: %.3e\n", max_err);

      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t8);
    }
  }

  if (do_bench) {
    struct timespec start, end;

    // Radix-2
    
    // Warm-up
    initialize_data_split(re, im, fftw_in_flat, n);
    fft_based_convolution(re, im, n, 2, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);
    memcpy(h_re_copy, h_re, n * sizeof(double));
    memcpy(h_im_copy, h_im, n * sizeof(double));
    permutation_avoiding_convolution(re, im, n, 2, twiddle_re, twiddle_im, h_re_copy, h_im_copy);
    
    // Timing runs
    double total_fft_based = 0., total_permutation_avoiding = 0.;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));

      clock_gettime(CLOCK_MONOTONIC, &start);
      fft_based_convolution(re, im, n, 2, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_fft_based += time_diff(start, end);

      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));

      clock_gettime(CLOCK_MONOTONIC, &start);
      permutation_avoiding_convolution(re, im, n, 2, twiddle_re, twiddle_im, h_re_copy, h_im_copy);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_permutation_avoiding += time_diff(start, end);
    }

    printf("Average radix-2 FFT-based time:               %.6f seconds\n", total_fft_based / timing_runs);
    printf("Average radix-2 permutation-avoiding time:    %.6f seconds\n", total_permutation_avoiding / timing_runs);
    
    // Clean-up
    free(rho);
    free_twiddles(twiddle_re, twiddle_im, t2);  

    // Radix-4
    if (supports_radix_4(t2)) {
      int t4 = t2 >> 1;
      int* rho = precompute_index_reversal_permutation_r4(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r4(n, &twiddle_re, &twiddle_im);

      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));

      // Warm-up
      initialize_data_split(re, im, fftw_in_flat, n);
      fft_based_convolution(re, im, n, 4, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);
      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));
      permutation_avoiding_convolution(re, im, n, 4, twiddle_re, twiddle_im, h_re_copy, h_im_copy);

      // Timing runs
      total_fft_based = 0., total_permutation_avoiding = 0.;
      for (int t = 0; t < timing_runs; ++t) {
        initialize_data_split(re, im, fftw_in_flat, n);

        memcpy(h_re_copy, h_re, n * sizeof(double));
        memcpy(h_im_copy, h_im, n * sizeof(double));

        clock_gettime(CLOCK_MONOTONIC, &start);
        fft_based_convolution(re, im, n, 4, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_fft_based += time_diff(start, end);

        memcpy(h_re_copy, h_re, n * sizeof(double));
        memcpy(h_im_copy, h_im, n * sizeof(double));

        clock_gettime(CLOCK_MONOTONIC, &start);
        permutation_avoiding_convolution(re, im, n, 4, twiddle_re, twiddle_im, h_re_copy, h_im_copy);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_permutation_avoiding += time_diff(start, end);
      }

      printf("Average radix-4 FFT-based time:               %.6f seconds\n", total_fft_based / timing_runs);
      printf("Average radix-4 permutation-avoiding time:    %.6f seconds\n", total_permutation_avoiding / timing_runs);

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

      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));
          
      // Warm-up
      initialize_data_split(re, im, fftw_in_flat, n);
      fft_based_convolution(re, im, n, 8, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);
      memcpy(h_re_copy, h_re, n * sizeof(double));
      memcpy(h_im_copy, h_im, n * sizeof(double));
      permutation_avoiding_convolution(re, im, n, 8, twiddle_re, twiddle_im, h_re_copy, h_im_copy);

      // Timing runs
      total_fft_based = 0., total_permutation_avoiding = 0.;
      for (int t = 0; t < timing_runs; ++t) {
        initialize_data_split(re, im, fftw_in_flat, n);

        memcpy(h_re_copy, h_re, n * sizeof(double));
        memcpy(h_im_copy, h_im, n * sizeof(double));

        clock_gettime(CLOCK_MONOTONIC, &start);
        fft_based_convolution(re, im, n, 8, twiddle_re, twiddle_im, rho, h_re_copy, h_im_copy);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_fft_based += time_diff(start, end);

        memcpy(h_re_copy, h_re, n * sizeof(double));
        memcpy(h_im_copy, h_im, n * sizeof(double));

        clock_gettime(CLOCK_MONOTONIC, &start);
        permutation_avoiding_convolution(re, im, n, 8, twiddle_re, twiddle_im, h_re_copy, h_im_copy);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_permutation_avoiding += time_diff(start, end);
      }

      printf("Average radix-8 FFT-based time:               %.6f seconds\n", total_fft_based / timing_runs);
      printf("Average radix-8 permutation-avoiding time:    %.6f seconds\n", total_permutation_avoiding / timing_runs);
    
      // Clean-up
      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t8); 
    }

    // FFTW

    // Warm-up
    initialize_data_split(re, im, fftw_in_flat, n);
    fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);

    // Timing runs
    double total_fftw = 0.;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fftw_based_convolution(p_forward, p_backward, fftw_in_flat, fftw_out_flat, h_re, h_im, tmp_re, tmp_im, n);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_fftw += time_diff(start, end);
    }
    
    printf("Average FFTW time:                            %.6f seconds\n", total_fftw / timing_runs);
  }

  // Clean-up
  fftw_destroy_plan(p_forward);
  fftw_destroy_plan(p_backward);
  fftw_free(fftw_in);
  fftw_free(fftw_out);
  free(re); free(im); 
  free(tmp_re); free(tmp_im); 
  free(h_re); free(h_im); 
  free(h_re_copy); free(h_im_copy); 

  return 0;
}