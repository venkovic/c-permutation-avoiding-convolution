/*
 * Copyright (c) 2025 Nicolas Venkovic
 * 
 * This file is part of c-permutation-avoiding-convolution.
 * 
 * This file is licensed under the MIT License.
 * For the full license text, see the LICENSE file in the root directory of this project.
 */

#include "fft.h"

void fft_2d(double* restrict re, double* restrict im, 
            int n1, int n2, int r, 
            double** twiddle_re1, double** twiddle_im1, 
            double** twiddle_re2, double** twiddle_im2, 
            const int* rho) {
  int n = n1 * n2;
  apply_permutation_split(re, im, n, rho);

  if (r == 2) {
    apply_butterflies_r2_2d(re, im, n1, n2,
                            (const double**) twiddle_re1, (const double**) twiddle_im1,
                            (const double**) twiddle_re2, (const double**) twiddle_im2);
  }
  else if (r == 4) {
    apply_butterflies_r4_2d(re, im, n1, n2,
                            (const double**) twiddle_re1, (const double**) twiddle_im1,
                            (const double**) twiddle_re2, (const double**) twiddle_im2);
  }  
}

void fft2_2d(double* restrict re, double* restrict im, 
             int n1, int n2, int r, 
             double** twiddle_re1, double** twiddle_im1, 
             double** twiddle_re2, double** twiddle_im2, 
             const int* rho) {
  int n = n1 * n2;

  if (r == 2) {
    apply_transposed_butterflies_r2_2d(re, im, n1, n2,
                                       (const double**) twiddle_re1, (const double**) twiddle_im1,
                                       (const double**) twiddle_re2, (const double**) twiddle_im2);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4_2d(re, im, n1, n2,
                                       (const double**) twiddle_re1, (const double**) twiddle_im1,
                                       (const double**) twiddle_re2, (const double**) twiddle_im2);
  }  

  apply_permutation_split(re, im, n, rho);
}

void ufft_2d(double* restrict re, double* restrict im, 
             int n1, int n2, int r, 
             double** twiddle_re1, double** twiddle_im1, 
             double** twiddle_re2, double** twiddle_im2) {
  if (r == 2) {
    apply_butterflies_r2_2d(re, im, n1, n2,
                            (const double**) twiddle_re1, (const double**) twiddle_im1,
                            (const double**) twiddle_re2, (const double**) twiddle_im2);

   }
  else if (r == 4) {
    apply_butterflies_r4_2d(re, im, n1, n2,
                            (const double**) twiddle_re1, (const double**) twiddle_im1,
                            (const double**) twiddle_re2, (const double**) twiddle_im2);

   }
}

void ifft_2d(double* restrict re, double* restrict im, 
             int n1, int n2, int r, 
             double** twiddle_re1, double** twiddle_im1, 
             double** twiddle_re2, double** twiddle_im2, 
             const int* rho) {
  int n = n1 * n2;

  apply_permutation_split(re, im, n, rho);

  if (r == 2) {
    apply_conjugate_butterflies_r2_2d(re, im, n1, n2,
                                      (const double**) twiddle_re1, (const double**) twiddle_im1,
                                      (const double**) twiddle_re2, (const double**) twiddle_im2);
  }
  else if (r == 4) {
    apply_conjugate_butterflies_r4_2d(re, im, n1, n2,
                                      (const double**) twiddle_re1, (const double**) twiddle_im1,
                                      (const double**) twiddle_re2, (const double**) twiddle_im2);
  }  

  for (int i = 0; i < n; i++)
    re[i] /= n, im[i] /= n;
}

void uifft_2d(double* restrict re, double* restrict im, 
              int n1, int n2, int r, 
              double** twiddle_re1, double** twiddle_im1, 
              double** twiddle_re2, double** twiddle_im2) {
  int n = n1 * n2;

  if (r == 2) {
    apply_conjugate_butterflies_r2_2d(re, im, n1, n2,
                                      (const double**) twiddle_re1, (const double**) twiddle_im1,
                                      (const double**) twiddle_re2, (const double**) twiddle_im2);
  }
  else if (r == 4) {
    apply_conjugate_butterflies_r4_2d(re, im, n1, n2,
                                      (const double**) twiddle_re1, (const double**) twiddle_im1,
                                      (const double**) twiddle_re2, (const double**) twiddle_im2);
  }  

  for (int i = 0; i < n; i++)
    re[i] /= n, im[i] /= n;
}