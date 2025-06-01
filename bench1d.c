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
  double complex* x  = malloc(n * sizeof(double complex));
  double* re  = malloc(n * sizeof(double));
  double* im  = malloc(n * sizeof(double));
  fftw_complex* fftw_in  = fftw_malloc(n * sizeof(fftw_complex));
  fftw_complex* fftw_out = fftw_malloc(n * sizeof(fftw_complex));
  double* fftw_in_flat  = (double*)fftw_in;
  double* fftw_out_flat = (double*)fftw_out;

  int* rho = precompute_index_reversal_permutation_r2(n);
  double** twiddle_re, **twiddle_im;
  precompute_twiddles_r2(n, &twiddle_re, &twiddle_im);
  fftw_plan p = fftw_plan_dft_1d(n, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);

  if (do_verify) {
    initialize_data(x, fftw_in_flat, n);

    fft(x, re, im, n, 2, twiddle_re, twiddle_im, rho);
    fftw_execute(p);

    double max_err = get_max_error(re, im, fftw_out_flat, n);
    printf("Max error FFT radix-2 vs FFTW: %.3e\n", max_err);

    free(rho);
    free_twiddles(twiddle_re, twiddle_im, t2);

    if (n % 4 == 0) {
      int t4 = t2 >> 1;
      int* rho = precompute_index_reversal_permutation_r4(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r4(n, &twiddle_re, &twiddle_im);
          
      initialize_data(x, fftw_in_flat, n);
        
      fft(x, re, im, n, 4, twiddle_re, twiddle_im, rho);
      fftw_execute(p);

      max_err = get_max_error(re, im, fftw_out_flat, n);
      printf("Max error FFT radix-4 vs FFTW: %.3e\n", max_err);          

      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t4);
    }

    if (n % 8 == 0) {
      int t8 = t2 >> 2;
      int* rho = precompute_index_reversal_permutation_r8(n);
      double** twiddle_re, **twiddle_im;
      precompute_twiddles_r8(n, &twiddle_re, &twiddle_im);
          
      initialize_data(x, fftw_in_flat, n);
        
      fft(x, re, im, n, 8, twiddle_re, twiddle_im, rho);
      fftw_execute(p);

      max_err = get_max_error(re, im, fftw_out_flat, n);
      printf("Max error FFT radix-8 vs FFTW: %.3e\n", max_err);

      free(rho);
      free_twiddles(twiddle_re, twiddle_im, t8);      
    }
  }

  if (do_bench) {
    struct timespec start, end;

    // Warm-up
    initialize_data(x, fftw_in_flat, n);
    fft(x, re, im, n, 2, twiddle_re, twiddle_im, rho);
    fftw_execute(p);

    // Timed FFT
    double total_fft = 0.0, total_perm = 0.0, total_ufft = 0.0, total_fftw = 0.0;
    for (int t = 0; t < timing_runs; ++t) {
      initialize_data(x, fftw_in_flat, n);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fft(x, re, im, n, 2, twiddle_re, twiddle_im, rho);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_fft += time_diff(start, end);

      clock_gettime(CLOCK_MONOTONIC, &start);
      apply_permutation(x, n, rho);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_perm += time_diff(start, end);

      clock_gettime(CLOCK_MONOTONIC, &start);
      ufft(re, im, n, 2, twiddle_re, twiddle_im);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_ufft += time_diff(start, end);

      clock_gettime(CLOCK_MONOTONIC, &start);
      fftw_execute(p);
      clock_gettime(CLOCK_MONOTONIC, &end);
      total_fftw += time_diff(start, end);
    }

    printf("Average FFT time:            %.6f seconds\n", total_fft / timing_runs);
    printf("Average permutation time:    %.6f seconds\n", total_perm / timing_runs);
    printf("Average Unordered FFT time:  %.6f seconds\n", total_ufft / timing_runs);
    printf("Average FFTW time:           %.6f seconds\n", total_fftw / timing_runs);

    free(rho);
    free_twiddles(twiddle_re, twiddle_im, t2);
  }

  // Cleanup
  fftw_destroy_plan(p);
  fftw_free(fftw_in);
  fftw_free(fftw_out);
  free(x); free(re); free(im); 

  return 0;
}

















