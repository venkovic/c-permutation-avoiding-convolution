#include "fft.h"

void fft(double complex* x, double* restrict re, double* restrict im, int n, 
         int r, double** twiddle_re, double** twiddle_im, const int* rho) {
  apply_permutation(x, n, rho);

  for (int i = 0; i < n; ++i) {
    double complex z = ((double complex*)x)[i];
    re[i] = creal(z);
    im[i] = cimag(z);
  }

  if (r == 2) {
    apply_butterflies_r2(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 4) {
    apply_butterflies_r4(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
  else if (r == 8) {
    //apply_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
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
    //apply_butterflies_r8(re, im, n, (const double**)twiddle_re, (const double**)twiddle_im);
  }
}
