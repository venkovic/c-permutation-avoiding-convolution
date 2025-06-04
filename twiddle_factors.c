/*
 * Copyright (c) 2025 Nicolas Venkovic
 * 
 * This file is part of c-permutation-avoiding-convolution.
 * 
 * This file is licensed under the MIT License.
 * For the full license text, see the LICENSE file in the root directory of this project.
 */

#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void precompute_twiddles_r2(int n, double*** twiddle_re, double*** twiddle_im) {
  int t = __builtin_ctz(n); // t = log2 n

  *twiddle_re = malloc(t * sizeof(double*));
  *twiddle_im = malloc(t * sizeof(double*));

  for (int q = 0; q < t; ++q) {
    int k = 1 << (q + 1);
    int ell = k / 2;

    (*twiddle_re)[q] = malloc(ell * sizeof(double));
    (*twiddle_im)[q] = malloc(ell * sizeof(double));

    for (int j = 0; j < ell; ++j) {
      // Twiddle factor: w_k^j
      double angle = -2.0 * M_PI * j / k;
      (*twiddle_re)[q][j] = cos(angle);
      (*twiddle_im)[q][j] = sin(angle);
    }
  }
}

void precompute_twiddles_r4(int n, double*** twiddle_re, double*** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 2;  // Divide by 4
    t++;
  } // t = log4 n

  *twiddle_re = malloc(t * sizeof(double*));
  *twiddle_im = malloc(t * sizeof(double*));

  for (int q = 0; q < t; ++q) {
    int k = 1 << (2 * (q + 1));  // k = 4^(q+1)
    int ell = k / 4;          

    (*twiddle_re)[q] = malloc(3 * ell * sizeof(double));
    (*twiddle_im)[q] = malloc(3 * ell * sizeof(double));

    for (int j = 0; j < ell; ++j) {
      // First twiddle factor: w_k^j
      double angle1 = -2.0 * M_PI * j / k;
      (*twiddle_re)[q][j] = cos(angle1);
      (*twiddle_im)[q][j] = sin(angle1);

      // Second twiddle factor: w_k^(2j)
      double angle2 = -2.0 * M_PI * (2 * j) / k;
      (*twiddle_re)[q][j + ell] = cos(angle2);
      (*twiddle_im)[q][j + ell] = sin(angle2);

      // Third twiddle factor: w_k^(3j)
      double angle3 = -2.0 * M_PI * (3 * j) / k;
      (*twiddle_re)[q][j + 2 * ell] = cos(angle3);
      (*twiddle_im)[q][j + 2 * ell] = sin(angle3);
    }
  }
}

void precompute_twiddles_r8(int n, double*** twiddle_re, double*** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 3;  // Divide by 8
    t++;
  } // t = log8 n

  *twiddle_re = malloc(t * sizeof(double*));
  *twiddle_im = malloc(t * sizeof(double*));

  for (int q = 0; q < t; ++q) {
    int k = 1 << (3 * (q + 1));  // k = 8^(q+1)
    int ell = k / 8;          

    (*twiddle_re)[q] = malloc(7 * ell * sizeof(double));
    (*twiddle_im)[q] = malloc(7 * ell * sizeof(double));

    for (int j = 0; j < ell; ++j) {
      // First twiddle factor: w_k^j
      double angle1 = -2.0 * M_PI * j / k;
      (*twiddle_re)[q][j] = cos(angle1);
      (*twiddle_im)[q][j] = sin(angle1);

      // Second twiddle factor: w_k^(2j)
      double angle2 = -2.0 * M_PI * (2 * j) / k;
      (*twiddle_re)[q][j + ell] = cos(angle2);
      (*twiddle_im)[q][j + ell] = sin(angle2);

      // Third twiddle factor: w_k^(3j)
      double angle3 = -2.0 * M_PI * (3 * j) / k;
      (*twiddle_re)[q][j + 2 * ell] = cos(angle3);
      (*twiddle_im)[q][j + 2 * ell] = sin(angle3);

      // Third twiddle factor: w_k^(4j)
      double angle4 = -2.0 * M_PI * (4 * j) / k;
      (*twiddle_re)[q][j + 3 * ell] = cos(angle4);
      (*twiddle_im)[q][j + 3 * ell] = sin(angle4);

      // Third twiddle factor: w_k^(5j)
      double angle5 = -2.0 * M_PI * (5 * j) / k;
      (*twiddle_re)[q][j + 4 * ell] = cos(angle5);
      (*twiddle_im)[q][j + 4 * ell] = sin(angle5);

      // Third twiddle factor: w_k^(6j)
      double angle6 = -2.0 * M_PI * (6 * j) / k;
      (*twiddle_re)[q][j + 5 * ell] = cos(angle6);
      (*twiddle_im)[q][j + 5 * ell] = sin(angle6);

      // Third twiddle factor: w_k^(7j)
      double angle7 = -2.0 * M_PI * (7 * j) / k;
      (*twiddle_re)[q][j + 6 * ell] = cos(angle7);
      (*twiddle_im)[q][j + 6 * ell] = sin(angle7);
    }
  }
}

void free_twiddles(double** twiddle_re, double** twiddle_im, int t) {
  for (int q = 0; q < t; ++q) {
    free(twiddle_re[q]);
    free(twiddle_im[q]);
  }
  free(twiddle_re);
  free(twiddle_im);
}