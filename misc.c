#include <time.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>

double time_diff(struct timespec start, struct timespec end) {
  return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

double get_max_error(double *re, double *im, double *fftw_flat, int n) {
  double max_err = 0.0;
  for (int i = 0; i < n; ++i) {
    double err_re = re[i] - fftw_flat[2 * i];
    double err_im = im[i] - fftw_flat[2 * i + 1];
    double err = hypot(err_re, err_im);
    if (err > max_err)
      max_err = err;
  }
  return max_err;
}

void initialize_data(double complex *x, double *fftw_in_flat, int n) {
  for (int i = 0; i < n; ++i) {
    x[i] = rand() / (double)RAND_MAX + I * rand() / (double)RAND_MAX;
    fftw_in_flat[2 * i]     = creal(x[i]);
    fftw_in_flat[2 * i + 1] = cimag(x[i]);
  }
}

void initialize_data_split(double *re, double *im, double *fftw_in_flat, int n) {
  for (int i = 0; i < n; ++i) {
    re[i] = rand() / (double)RAND_MAX;
    im[i] = rand() / (double)RAND_MAX; 
    fftw_in_flat[2 * i]     = re[i];
    fftw_in_flat[2 * i + 1] = im[i];
  }
}

void initialize_filter(double *g_re, double *g_im, int n) {
  for (int i = 0; i < n; ++i) {
    g_re[i] = rand() / (double)RAND_MAX;
    g_im[i] = rand() / (double)RAND_MAX; 
  }
}

bool supports_radix_4(int t2) {
  if (t2 % 2) return false;
  return true;
}

bool supports_radix_8(int t2) {
  if (t2 % 3) return false;
  return true;
}