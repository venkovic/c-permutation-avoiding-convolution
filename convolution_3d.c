/*
 * Copyright (c) 2025 Nicolas Venkovic
 * 
 * This file is part of c-permutation-avoiding-convolution.
 * 
 * This file is licensed under the MIT License.
 * For the full license text, see the LICENSE file in the root directory of this project.
 */

#include "fft.h"

void fft_based_convolution_3d(double* restrict re, double* restrict im, 
                              int n1, int n2, int n3, int r, 
                              double** twiddle_re1, double** twiddle_im1, 
                              double** twiddle_re2, double** twiddle_im2,
                              double** twiddle_re3, double** twiddle_im3,                               
                              const int* rho, 
                              double* restrict g_re, double* restrict g_im) {
  int n = n1 * n2 * n3;

  fft_3d(re, im, n1, n2, n3, r, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, twiddle_re3, twiddle_im3, rho);

  for (int i = 0; i < n; i++) {
    double x_re = re[i], x_im = im[i];
    re[i] = x_re * g_re[i] - x_im * g_im[i];
    im[i] = x_re * g_im[i] + x_im * g_re[i];
  }

  ifft_3d(re, im, n1, n2, n3, r, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, twiddle_re3, twiddle_im3, rho);
}

void permutation_avoiding_convolution_3d(double* restrict re, double* restrict im, 
                                         int n1, int n2, int n3, int r, 
                                         double** twiddle_re1, double** twiddle_im1,
                                         double** twiddle_re2, double** twiddle_im2,
                                         double** twiddle_re3, double** twiddle_im3,
                                         double* restrict Pg_re, double* restrict Pg_im) {
  int n = n1 * n2 * n3;

  if (r == 2) {
    apply_transposed_butterflies_r2_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }

  for (int i = 0; i < n; i++) {
    double x_re = re[i], x_im = im[i];
    re[i] = (x_re * Pg_re[i] - x_im * Pg_im[i]) / n;
    im[i] = (x_re * Pg_im[i] + x_im * Pg_re[i]) / n;
  }

  if (r == 2) {
    apply_conjugate_butterflies_r2_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }
  else if (r == 4) {
    apply_conjugate_butterflies_r4_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }
}