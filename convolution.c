/*
 * Copyright (c) 2025 Nicolas Venkovic
 * 
 * This file is part of c-permutation-avoiding-convolution.
 * 
 * This file is licensed under the MIT License.
 * For the full license text, see the LICENSE file in the root directory of this project.
 */

#include "fft.h"

void fft_based_convolution(double* restrict re, double* restrict im, int n, int r, 
                           double** twiddle_re, double** twiddle_im, const int* rho, 
                           double* restrict h_re, double* restrict h_im) {

  fft(re, im, n, r, twiddle_re, twiddle_im, rho);

  fft(h_re, h_im, n, r, twiddle_re, twiddle_im, rho);

  for (int i = 0; i < n; i++) {
    double Fx_re = re[i], Fx_im = im[i];
    double Fh_re = h_re[i], Fh_im = h_im[i];

    re[i] = Fx_re * Fh_re - Fx_im * Fh_im;
    im[i] = Fx_re * Fh_im + Fx_im * Fh_re;
  }

  ifft(re, im, n, r, twiddle_re, twiddle_im, rho);
}

void permutation_avoiding_convolution(double* restrict re, double* restrict im, int n, 
                                      int r, double** twiddle_re, double** twiddle_im, 
                                      double* restrict h_re, double* restrict h_im) {

  if (r == 2) {
    apply_transposed_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_transposed_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }

  if (r == 2) {
    apply_transposed_butterflies_r2(h_re, h_im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4(h_re, h_im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_transposed_butterflies_r8(h_re, h_im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }

  for (int i = 0; i < n; i++) {
    double Atx_re = re[i], Atx_im = im[i];
    double Ath_re = h_re[i], Ath_im = h_im[i];

    re[i] = (Atx_re * Ath_re - Atx_im * Ath_im) / n;
    im[i] = (Atx_re * Ath_im + Atx_im * Ath_re) / n;
  }

  if (r == 2) {
    apply_conjugate_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_conjugate_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_conjugate_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
}