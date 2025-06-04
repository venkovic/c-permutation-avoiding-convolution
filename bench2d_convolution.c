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
                            double* fftw_out_flat,
                            double* g_re, double* g_im, 
                            int n) {
  
  fftw_execute(p_forward);

  for (int i = 0; i < n; i++) {
    double re = fftw_out_flat[2 * i];
    double im = fftw_out_flat[2 * i + 1];

    fftw_out_flat[2 * i]     = (re * g_re[i] - im * g_im[i]) / n;
    fftw_out_flat[2 * i + 1] = (re * g_im[i] + im * g_re[i]) / n;
  }

  fftw_execute(p_backward);
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
  double* g_re  = malloc(n * sizeof(double));
  double* g_im  = malloc(n * sizeof(double));
  double* Pg_re  = malloc(n * sizeof(double));
  double* Pg_im  = malloc(n * sizeof(double));
  fftw_complex* fftw_in  = fftw_malloc(n * sizeof(fftw_complex));
  fftw_complex* fftw_out = fftw_malloc(n * sizeof(fftw_complex));
  double* fftw_in_flat  = (double*)fftw_in;
  double* fftw_out_flat = (double*)fftw_out;

  int* rho = precompute_index_reversal_permutation_r2_2d(n1, n2);
  double** twiddle_re1, **twiddle_im1;
  double** twiddle_re2, **twiddle_im2;
  precompute_twiddles_r2(n1, &twiddle_re1, &twiddle_im1);
  precompute_twiddles_r2(n2, &twiddle_re2, &twiddle_im2);

  fftw_plan p_forward = fftw_plan_dft_2d(n1, n2, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan p_backward = fftw_plan_dft_2d(n1, n2, fftw_out, fftw_in, FFTW_BACKWARD, FFTW_ESTIMATE);

  initialize_filter(g_re, g_im, n);
  memcpy(Pg_re, g_re, n * sizeof(double));
  memcpy(Pg_im, g_im, n * sizeof(double));
  apply_permutation_split(Pg_re, Pg_im, n, rho);

  if (do_verify) {
    initialize_data_split(re, im, fftw_in_flat, n);
    
    // FFTW
    fftw_based_convolution(p_forward, p_backward, fftw_out_flat, g_re, g_im, n);

    // Radix-2
    fft_based_convolution_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho, g_re, g_im);

    double max_err = get_max_error(re, im, fftw_in_flat, n);
    printf("Max error FFT-based radix-2 vs FFTW: %.3e\n", max_err);

    initialize_data_split(re, im, fftw_in_flat, n);
    
    // FFTW
    fftw_based_convolution(p_forward, p_backward, fftw_out_flat, g_re, g_im, n);

    // Permutation-avoiding radix-2
    permutation_avoiding_convolution_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, Pg_re, Pg_im);

    max_err = get_max_error(re, im, fftw_in_flat, n);
    printf("Max error Permutation-avoiding radix-2 vs FFTW: %.3e\n", max_err);

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

      memcpy(Pg_re, g_re, n * sizeof(double));
      memcpy(Pg_im, g_im, n * sizeof(double));
      apply_permutation_split(Pg_re, Pg_im, n, rho);

      initialize_data_split(re, im, fftw_in_flat, n);

      // FFTW
      fftw_based_convolution(p_forward, p_backward, fftw_out_flat, g_re, g_im, n);

      // FFT-based Radix-4
      fft_based_convolution_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho, g_re, g_im);

      max_err = get_max_error(re, im, fftw_in_flat, n);
      printf("Max error FFT-based radix-4 vs FFTW: %.3e\n", max_err);      
       
      initialize_data_split(re, im, fftw_in_flat, n);   

      // FFTW
      fftw_based_convolution(p_forward, p_backward, fftw_out_flat, g_re, g_im, n);

      // Permutation-avoiding radix-4
      permutation_avoiding_convolution_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, Pg_re, Pg_im);

      max_err = get_max_error(re, im, fftw_in_flat, n);
      printf("Max error Permutation-avoiding radix-4 vs FFTW: %.3e\n", max_err);

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
    fft_based_convolution_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho, g_re, g_im);
    permutation_avoiding_convolution_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, Pg_re, Pg_im);
    
    // Timing runs
    double total_fft_based = 0., total_permutation_avoiding = 0.;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fft_based_convolution_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho, g_re, g_im);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_fft_based += time_diff(start, end);

      clock_gettime(CLOCK_MONOTONIC, &start);
      permutation_avoiding_convolution_2d(re, im, n1, n2, 2, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, Pg_re, Pg_im);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_permutation_avoiding += time_diff(start, end);
    }

    printf("Average radix-2 FFT-based time:               %.6f seconds\n", total_fft_based / timing_runs);
    printf("Average radix-2 permutation-avoiding time:    %.6f seconds\n", total_permutation_avoiding / timing_runs);

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

      memcpy(Pg_re, g_re, n * sizeof(double));
      memcpy(Pg_im, g_im, n * sizeof(double));
      apply_permutation_split(Pg_re, Pg_im, n, rho);

      // Warm-up
      initialize_data_split(re, im, fftw_in_flat, n);
      fft_based_convolution_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho, g_re, g_im);
      permutation_avoiding_convolution_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, Pg_re, Pg_im);
        
      // Timing runs
      total_fft_based = 0., total_permutation_avoiding = 0.;
      for (int t = 0; t < timing_runs; ++t) {
        initialize_data_split(re, im, fftw_in_flat, n);

        clock_gettime(CLOCK_MONOTONIC, &start);
        fft_based_convolution_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, rho, g_re, g_im);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_fft_based += time_diff(start, end);

        clock_gettime(CLOCK_MONOTONIC, &start);
        permutation_avoiding_convolution_2d(re, im, n1, n2, 4, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, Pg_re, Pg_im);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_permutation_avoiding += time_diff(start, end);
      }

      printf("Average radix-4 FFT-based time:               %.6f seconds\n", total_fft_based / timing_runs);
      printf("Average radix-4 permutation-avoiding time:    %.6f seconds\n", total_permutation_avoiding / timing_runs);

      // Clean-up
      free(rho);
      free_twiddles(twiddle_re1, twiddle_im1, t4_1); 
      free_twiddles(twiddle_re2, twiddle_im2, t4_2); 
    }

    // FFTW

    // Warm-up
    initialize_data_split(re, im, fftw_in_flat, n);
    fftw_based_convolution(p_forward, p_backward, fftw_out_flat, g_re, g_im, n);

    // Timing runs
    double total_fftw = 0.;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data_split(re, im, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fftw_based_convolution(p_forward, p_backward, fftw_out_flat, g_re, g_im, n);
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

  return 0;
}