#ifndef FFT_H
#define FFT_H

#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>

// In twiddle_factors.c
void precompute_twiddles_r2(int n, double*** twiddle_re, double*** twiddle_im);
void precompute_twiddles_r4(int n, double*** twiddle_re, double*** twiddle_im);
void precompute_twiddles_r8(int n, double*** twiddle_re, double*** twiddle_im);
void free_twiddles(double** twiddle_re, double** twiddle_im, int t);

// In butterflies.c
void apply_butterflies_r2(double* restrict re, double* restrict im, int n,
                          const double** twiddle_re, const double** twiddle_im);
void apply_conjugate_butterflies_r2(double* restrict re, double* restrict im, int n,
                                    const double** twiddle_re, const double** twiddle_im);
void apply_transposed_butterflies_r2(double* restrict re, double* restrict im, int n,
                                     const double** twiddle_re, const double** twiddle_im);
void apply_butterflies_r4(double* restrict re, double* restrict im, int n,
                          const double** twiddle_re, const double** twiddle_im);
void apply_conjugate_butterflies_r4(double* restrict re, double* restrict im, int n,
                                    const double** twiddle_re, const double** twiddle_im);
void apply_transposed_butterflies_r4(double* restrict re, double* restrict im, int n,
                                     const double** twiddle_re, const double** twiddle_im);
void apply_butterflies_r8(double* restrict re, double* restrict im, int n,
                          const double** twiddle_re, const double** twiddle_im);
void apply_conjugate_butterflies_r8(double* restrict re, double* restrict im, int n,
                                    const double** twiddle_re, const double** twiddle_im);
void apply_transposed_butterflies_r8(double* restrict re, double* restrict im, int n,
                                     const double** twiddle_re, const double** twiddle_im);

// In index_reversals.c
int reverse_index_r2(int x, int t);
int* precompute_index_reversal_permutation_r2(int n);
int reverse_index_r4(int x, int t);
int* precompute_index_reversal_permutation_r4(int n);
int reverse_index_r8(int x, int t);
int* precompute_index_reversal_permutation_r8(int n);
void apply_permutation(double complex* x, int n, const int* rho);
void apply_permutation_split(double* restrict re, double* restrict im, 
                             int n, const int* rho);

// In fft.c
void fft(double* restrict re, double* restrict im, int n, int r, 
         double** twiddle_re, double** twiddle_im, const int* rho);
void fft2(double* restrict re, double* restrict im, int n, int r, 
          double** twiddle_re, double** twiddle_im, const int* rho);
void ufft(double* restrict re, double* restrict im, int n, 
          int r, double** twiddle_re, double** twiddle_im);
void ifft(double* restrict re, double* restrict im, int n, int r, 
          double** twiddle_re, double** twiddle_im, const int* rho);
void uifft(double* restrict re, double* restrict im, int n, 
           int r, double** twiddle_re, double** twiddle_im);

// In convolution.c
void fft_based_convolution(double* restrict re, double* restrict im, int n, int r, 
                           double** twiddle_re, double** twiddle_im, const int* rho, 
                           double* restrict g_re, double* restrict g_im);
void permutation_avoiding_convolution(double* restrict re, double* restrict im, int n, 
                                      int r, double** twiddle_re, double** twiddle_im, 
                                      double* restrict Pg_re, double* restrict Pg_im);

// In misc.c
double time_diff(struct timespec start, struct timespec end);
double get_max_error(double *re, double *im, double *fftw_flat, int n);
void initialize_data(double complex *x, double *fftw_in_flat, int n);
void initialize_data_split(double *re, double *im, double *fftw_in_flat, int n);
void initialize_filter(double *g_re, double *g_im, int n);
bool supports_radix_4(int t2);
bool supports_radix_8(int t2);

#endif // FFT_H