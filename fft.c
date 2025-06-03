#include "fft.h"

void fft(double* restrict re, double* restrict im, int n, int r, 
         double** twiddle_re, double** twiddle_im, const int* rho) {
  apply_permutation_split(re, im, n, rho);

  if (r == 2) {
    apply_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
}

void fft2(double* restrict re, double* restrict im, int n, int r, 
          double** twiddle_re, double** twiddle_im, const int* rho) {
  if (r == 2) {
    apply_transposed_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_transposed_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_transposed_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }

  apply_permutation_split(re, im, n, rho);
}

void ufft(double* restrict re, double* restrict im, int n, 
          int r, double** twiddle_re, double** twiddle_im) {
  if (r == 2) {
    apply_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
}

void ifft(double* restrict re, double* restrict im, int n, int r, 
          double** twiddle_re, double** twiddle_im, const int* rho) {
  apply_permutation_split(re, im, n, rho);

  if (r == 2) {
    apply_conjugate_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_conjugate_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_conjugate_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }

  for (int i = 0; i < n; i++)
    re[i] /= n, im[i] /= n;
}

void uifft(double* restrict re, double* restrict im, int n, 
           int r, double** twiddle_re, double** twiddle_im) {

  if (r == 2) {
    apply_conjugate_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_conjugate_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    apply_conjugate_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }

  for (int i = 0; i < n; i++)
    re[i] /= n, im[i] /= n;
}