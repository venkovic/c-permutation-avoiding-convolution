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
                              double* restrict h_re, double* restrict h_im) {
  int n = n1 * n2 * n3;

  fft_3d(re, im, n1, n2, n3, r, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, twiddle_re3, twiddle_im3, rho);

  fft_3d(h_re, h_im, n1, n2, n3, r, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, twiddle_re3, twiddle_im3, rho);

  for (int i = 0; i < n; i++) {
    double Fx_re = re[i], Fx_im = im[i];
    double Fh_re = h_re[i], Fh_im = h_im[i];

    re[i] = Fx_re * Fh_re - Fx_im * Fh_im;
    im[i] = Fx_re * Fh_im + Fx_im * Fh_re;
  }

  ifft_3d(re, im, n1, n2, n3, r, twiddle_re1, twiddle_im1, twiddle_re2, twiddle_im2, twiddle_re3, twiddle_im3, rho);
}

void permutation_avoiding_convolution_3d(double* restrict re, double* restrict im, 
                                         int n1, int n2, int n3, int r, 
                                         double** twiddle_re1, double** twiddle_im1,
                                         double** twiddle_re2, double** twiddle_im2,
                                         double** twiddle_re3, double** twiddle_im3,
                                         double* restrict h_re, double* restrict h_im) {
  int n = n1 * n2 * n3;

  if (r == 2) {
    apply_transposed_butterflies_r2_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }

  if (r == 2) {
    apply_transposed_butterflies_r2_3d(h_re, h_im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4_3d(h_re, h_im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }

  for (int i = 0; i < n; i++) {
    double Atx_re = re[i], Atx_im = im[i];
    double Ath_re = h_re[i], Ath_im = h_im[i];

    re[i] = (Atx_re * Ath_re - Atx_im * Ath_im) / n;
    im[i] = (Atx_re * Ath_im + Atx_im * Ath_re) / n;
  }

  if (r == 2) {
    apply_conjugate_butterflies_r2_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }
  else if (r == 4) {
    apply_conjugate_butterflies_r4_3d(re, im, n1, n2, n3, (const double**)twiddle_re1, (const double**)twiddle_im1, (const double**)twiddle_re2, (const double**)twiddle_im2, (const double**)twiddle_re3, (const double**)twiddle_im3);
  }
}