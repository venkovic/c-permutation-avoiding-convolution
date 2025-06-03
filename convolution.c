#include "fft.h"

void fft_based_convolution(double* restrict re, double* restrict im, int n, int r, 
                           double** twiddle_re, double** twiddle_im, const int* rho, 
                           double* restrict g_re, double* restrict g_im) {

  fft(re, im, n, r, twiddle_re, twiddle_im, rho);

  for (int i = 0; i < n; i++) {
    double x_re = re[i], x_im = im[i];
    re[i] = x_re * g_re[i] - x_im * g_im[i];
    im[i] = x_re * g_im[i] + x_im * g_re[i];
  }

  ifft(re, im, n, r, twiddle_re, twiddle_im, rho);
}

void permutation_avoiding_convolution(double* restrict re, double* restrict im, int n, 
                                      int r, double** twiddle_re, double** twiddle_im, 
                                      double* restrict Pg_re, double* restrict Pg_im) {

  if (r == 2) {
    apply_transposed_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_transposed_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }

  for (int i = 0; i < n; i++) {
    double x_re = re[i], x_im = im[i];
    re[i] = (x_re * Pg_re[i] - x_im * Pg_im[i]) / n;
    im[i] = (x_re * Pg_im[i] + x_im * Pg_re[i]) / n;
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