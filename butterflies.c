#include <math.h>

void apply_butterflies_r2(double* restrict re, double* restrict im, int n,
                          const double** twiddle_re, const double** twiddle_im) {
  int t = __builtin_ctz(n); // t = log2 n

  for (int q = 0; q < t; ++q) {
    int k = 1 << (q + 1);
    int ell = k / 2;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double tau_re = twiddle_re[q][j] * re[i + ell] - twiddle_im[q][j] * im[i + ell];
        double tau_im = twiddle_re[q][j] * im[i + ell] + twiddle_im[q][j] * re[i + ell];

        double z_re = re[i];
        double z_im = im[i];

        re[i]       = z_re + tau_re;
        im[i]       = z_im + tau_im;
        re[i + ell] = z_re - tau_re;
        im[i + ell] = z_im - tau_im;
      }
    }
  }
}

void apply_conjugate_butterflies_r2(double* restrict re, double* restrict im, int n,
                                    const double** twiddle_re, const double** twiddle_im) {
  int t = __builtin_ctz(n); // t = log2 n

  for (int q = 0; q < t; ++q) {
    int k = 1 << (q + 1);
    int ell = k / 2;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double tau_re = twiddle_re[q][j] * re[i + ell] + twiddle_im[q][j] * im[i + ell];
        double tau_im = twiddle_re[q][j] * im[i + ell] - twiddle_im[q][j] * re[i + ell];

        double z_re = re[i];
        double z_im = im[i];

        re[i]       = z_re + tau_re;
        im[i]       = z_im + tau_im;
        re[i + ell] = z_re - tau_re;
        im[i + ell] = z_im - tau_im;
      }
    }
  }
}

void apply_transposed_butterflies_r2(double* restrict re, double* restrict im, int n,
                                     const double** twiddle_re, const double** twiddle_im) {
  int t = __builtin_ctz(n); // t = log2 n

  for (int q = 0; q < t; ++q) {
    int k = 1 << (q + 1);
    int ell = k / 2;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double tau_re = re[i + ell];
        double tau_im = im[i + ell];

        re[i + ell] = twiddle_re[q][j] * (re[i] - tau_re) - twiddle_im[q][j] * (im[i] - tau_im);
        im[i + ell] = twiddle_re[q][j] * (im[i] - tau_im) + twiddle_im[q][j] * (re[i] - tau_re);
        re[i]       = re[i] + tau_re;
        im[i]       = im[i] + tau_im;
      }
    }
  }
}

void apply_butterflies_r4(double* restrict re, double* restrict im, int n,
                          const double** twiddle_re, const double** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 2; // Divide by 4
    t++;
  } // t = log4 n
    
  for (int q = 0; q < t; ++q) {
    int k = 1 << (2 * (q + 1));  // k = 4^(q+1)
    int ell = k / 4;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double x1_re = re[i];
        double x1_im = im[i];
        double x2_re = re[i + ell];
        double x2_im = im[i + ell];
        double x3_re = re[i + 2 * ell];
        double x3_im = im[i + 2 * ell];
        double x4_re = re[i + 3 * ell];
        double x4_im = im[i + 3 * ell];

        double z1_re = x1_re;
        double z1_im = x1_im;
        double z2_re = twiddle_re[q][j] * x2_re - twiddle_im[q][j] * x2_im;
        double z2_im = twiddle_re[q][j] * x2_im + twiddle_im[q][j] * x2_re;
        double z3_re = twiddle_re[q][j + ell] * x3_re - twiddle_im[q][j + ell] * x3_im;
        double z3_im = twiddle_re[q][j + ell] * x3_im + twiddle_im[q][j + ell] * x3_re;
        double z4_re = twiddle_re[q][j + 2 * ell] * x4_re - twiddle_im[q][j + 2 * ell] * x4_im;
        double z4_im = twiddle_re[q][j + 2 * ell] * x4_im + twiddle_im[q][j + 2 * ell] * x4_re;

        double tau1_re = z1_re + z3_re;
        double tau1_im = z1_im + z3_im;
        double tau2_re = z1_re - z3_re;
        double tau2_im = z1_im - z3_im;
        double tau3_re = z2_re + z4_re;
        double tau3_im = z2_im + z4_im;
        double tau4_re = z2_re - z4_re;
        double tau4_im = z2_im - z4_im;

        re[i]           = tau1_re + tau3_re;
        im[i]           = tau1_im + tau3_im;
        re[i + ell]     = tau2_re + tau4_im;
        im[i + ell]     = tau2_im - tau4_re;
        re[i + 2 * ell] = tau1_re - tau3_re;
        im[i + 2 * ell] = tau1_im - tau3_im;
        re[i + 3 * ell] = tau2_re - tau4_im;
        im[i + 3 * ell] = tau2_im + tau4_re;
      }
    }
  }
}

void apply_conjugate_butterflies_r4(double* restrict re, double* restrict im, int n,
                                    const double** twiddle_re, const double** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 2; // Divide by 4
    t++;
  } // t = log4 n
    
  for (int q = 0; q < t; ++q) {
    int k = 1 << (2 * (q + 1));  // k = 4^(q+1)
    int ell = k / 4;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double x1_re = re[i];
        double x1_im = im[i];
        double x2_re = re[i + ell];
        double x2_im = im[i + ell];
        double x3_re = re[i + 2 * ell];
        double x3_im = im[i + 2 * ell];
        double x4_re = re[i + 3 * ell];
        double x4_im = im[i + 3 * ell];

        double z1_re = x1_re;
        double z1_im = x1_im;
        double z2_re = twiddle_re[q][j] * x2_re + twiddle_im[q][j] * x2_im;
        double z2_im = twiddle_re[q][j] * x2_im - twiddle_im[q][j] * x2_re;
        double z3_re = twiddle_re[q][j + ell] * x3_re + twiddle_im[q][j + ell] * x3_im;
        double z3_im = twiddle_re[q][j + ell] * x3_im - twiddle_im[q][j + ell] * x3_re;
        double z4_re = twiddle_re[q][j + 2 * ell] * x4_re + twiddle_im[q][j + 2 * ell] * x4_im;
        double z4_im = twiddle_re[q][j + 2 * ell] * x4_im - twiddle_im[q][j + 2 * ell] * x4_re;

        double tau1_re = z1_re + z3_re;
        double tau1_im = z1_im + z3_im;
        double tau2_re = z1_re - z3_re;
        double tau2_im = z1_im - z3_im;
        double tau3_re = z2_re + z4_re;
        double tau3_im = z2_im + z4_im;
        double tau4_re = z2_re - z4_re;
        double tau4_im = z2_im - z4_im;

        re[i]           = tau1_re + tau3_re;
        im[i]           = tau1_im + tau3_im;
        re[i + ell]     = tau2_re - tau4_im;
        im[i + ell]     = tau2_im + tau4_re;
        re[i + 2 * ell] = tau1_re - tau3_re;
        im[i + 2 * ell] = tau1_im - tau3_im;
        re[i + 3 * ell] = tau2_re + tau4_im;
        im[i + 3 * ell] = tau2_im - tau4_re;
      }
    }
  }
}

void apply_transposed_butterflies_r4(double* restrict re, double* restrict im, int n,
                                     const double** twiddle_re, const double** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 2; // Divide by 4
    t++;
  } // t = log4 n
    
  for (int q = 0; q < t; ++q) {
    int k = 1 << (2 * (q + 1));  // k = 4^(q+1)
    int ell = k / 4;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double z1_re = re[i];
        double z1_im = im[i];
        double z2_re = re[i + ell];
        double z2_im = im[i + ell];
        double z3_re = re[i + 2 * ell];
        double z3_im = im[i + 2 * ell];
        double z4_re = re[i + 3 * ell];
        double z4_im = im[i + 3 * ell];

        double tau1_re = z1_re + z3_re;
        double tau1_im = z1_im + z3_im;
        double tau2_re = twiddle_re[q][j] * (z1_re - z3_re) - twiddle_im[q][j] * (z1_im - z3_im);
        double tau2_im = twiddle_re[q][j] * (z1_im - z3_im) + twiddle_im[q][j] * (z1_re - z3_re);
        double tau3_re = z2_re + z4_re;
        double tau3_im = z2_im + z4_im;  
        double tau4_re = twiddle_re[q][j] * (z2_re - z4_re) - twiddle_im[q][j] * (z2_im - z4_im);
        double tau4_im = twiddle_re[q][j] * (z2_im - z4_im) + twiddle_im[q][j] * (z2_re - z4_re);

        re[i]           = tau1_re + tau3_re;
        im[i]           = tau1_im + tau3_im;
        re[i + ell]     = tau2_re + tau4_im;
        im[i + ell]     = tau2_im - tau4_re;
        re[i + 2 * ell] = twiddle_re[q][j + ell] * (tau1_re - tau3_re) - twiddle_im[q][j + ell] * (tau1_im - tau3_im);
        im[i + 2 * ell] = twiddle_re[q][j + ell] * (tau1_im - tau3_im) + twiddle_im[q][j + ell] * (tau1_re - tau3_re);
        re[i + 3 * ell] = twiddle_re[q][j + ell] * (tau2_re - tau4_im) - twiddle_im[q][j + ell] * (tau2_im + tau4_re);
        im[i + 3 * ell] = twiddle_re[q][j + ell] * (tau2_im + tau4_re) + twiddle_im[q][j + ell] * (tau2_re - tau4_im);
      }
    }
  }
}

void apply_butterflies_r8(double* restrict re, double* restrict im, int n,
                          const double** twiddle_re, const double** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 3; // Divide by 8
    t++;
  } // t = log8 n

  double tmp = 1. / sqrt(2.);
  double a_re =  tmp;
  double a_im = -tmp;
  double b_re = -tmp;
  double b_im = -tmp;
    
  for (int q = 0; q < t; ++q) {
    int k = 1 << (3 * (q + 1));  // k = 8^(q+1)
    int ell = k / 8;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double x1_re = re[i];
        double x1_im = im[i];
        double x2_re = re[i + ell];
        double x2_im = im[i + ell];
        double x3_re = re[i + 2 * ell];
        double x3_im = im[i + 2 * ell];
        double x4_re = re[i + 3 * ell];
        double x4_im = im[i + 3 * ell];
        double x5_re = re[i + 4 * ell];
        double x5_im = im[i + 4 * ell];
        double x6_re = re[i + 5 * ell];
        double x6_im = im[i + 5 * ell];
        double x7_re = re[i + 6 * ell];
        double x7_im = im[i + 6 * ell];
        double x8_re = re[i + 7 * ell];
        double x8_im = im[i + 7 * ell];

        double z1_re = x1_re;
        double z1_im = x1_im;
        double z2_re = twiddle_re[q][j] * x2_re - twiddle_im[q][j] * x2_im;
        double z2_im = twiddle_re[q][j] * x2_im + twiddle_im[q][j] * x2_re;
        double z3_re = twiddle_re[q][j + ell] * x3_re - twiddle_im[q][j + ell] * x3_im;
        double z3_im = twiddle_re[q][j + ell] * x3_im + twiddle_im[q][j + ell] * x3_re;
        double z4_re = twiddle_re[q][j + 2 * ell] * x4_re - twiddle_im[q][j + 2 * ell] * x4_im;
        double z4_im = twiddle_re[q][j + 2 * ell] * x4_im + twiddle_im[q][j + 2 * ell] * x4_re;
        double z5_re = twiddle_re[q][j + 3 * ell] * x5_re - twiddle_im[q][j + 3 * ell] * x5_im;
        double z5_im = twiddle_re[q][j + 3 * ell] * x5_im + twiddle_im[q][j + 3 * ell] * x5_re;
        double z6_re = twiddle_re[q][j + 4 * ell] * x6_re - twiddle_im[q][j + 4 * ell] * x6_im;
        double z6_im = twiddle_re[q][j + 4 * ell] * x6_im + twiddle_im[q][j + 4 * ell] * x6_re;
        double z7_re = twiddle_re[q][j + 5 * ell] * x7_re - twiddle_im[q][j + 5 * ell] * x7_im;
        double z7_im = twiddle_re[q][j + 5 * ell] * x7_im + twiddle_im[q][j + 5 * ell] * x7_re;
        double z8_re = twiddle_re[q][j + 6 * ell] * x8_re - twiddle_im[q][j + 6 * ell] * x8_im;
        double z8_im = twiddle_re[q][j + 6 * ell] * x8_im + twiddle_im[q][j + 6 * ell] * x8_re;

        double tau1_re = z1_re + z5_re;
        double tau1_im = z1_im + z5_im;
        double tau2_re = z1_re - z5_re;
        double tau2_im = z1_im - z5_im;
        double tau3_re = z2_re + z6_re;
        double tau3_im = z2_im + z6_im;
        double tau4_re = z2_re - z6_re;
        double tau4_im = z2_im - z6_im;
        double tau5_re = z3_re + z7_re;
        double tau5_im = z3_im + z7_im;
        double tau6_re = z3_re - z7_re;
        double tau6_im = z3_im - z7_im;
        double tau7_re = z4_re + z8_re;
        double tau7_im = z4_im + z8_im;
        double tau8_re = z4_re - z8_re;
        double tau8_im = z4_im - z8_im;

        double a_re_x_tau4_re = a_re * tau4_re;
        double a_re_x_tau4_im = a_re * tau4_im;
        double a_im_x_tau4_re = a_im * tau4_re;
        double a_im_x_tau4_im = a_im * tau4_im;
        double a_re_x_tau8_re = a_re * tau8_re;
        double a_re_x_tau8_im = a_re * tau8_im;
        double a_im_x_tau8_re = a_im * tau8_re;
        double a_im_x_tau8_im = a_im * tau8_im;
        double b_re_x_tau4_re = b_re * tau4_re;
        double b_re_x_tau4_im = b_re * tau4_im;
        double b_im_x_tau4_re = b_im * tau4_re;
        double b_im_x_tau4_im = b_im * tau4_im;
        double b_re_x_tau8_re = b_re * tau8_re;
        double b_re_x_tau8_im = b_re * tau8_im;
        double b_im_x_tau8_re = b_im * tau8_re;
        double b_im_x_tau8_im = b_im * tau8_im;

        re[i]           = tau1_re + tau3_re + tau5_re + tau7_re;
        im[i]           = tau1_im + tau3_im + tau5_im + tau7_im;
        re[i + ell]     = tau2_re + a_re_x_tau4_re - a_im_x_tau4_im + tau6_im + b_re_x_tau8_re - b_im_x_tau8_im;
        im[i + ell]     = tau2_im + a_re_x_tau4_im + a_im_x_tau4_re - tau6_re + b_re_x_tau8_im + b_im_x_tau8_re;
        re[i + 2 * ell] = tau1_re + tau3_im - tau5_re - tau7_im;
        im[i + 2 * ell] = tau1_im - tau3_re - tau5_im + tau7_re;
        re[i + 3 * ell] = tau2_re + b_re_x_tau4_re - b_im_x_tau4_im - tau6_im + a_re_x_tau8_re - a_im_x_tau8_im;
        im[i + 3 * ell] = tau2_im + b_re_x_tau4_im + b_im_x_tau4_re + tau6_re + a_re_x_tau8_im + a_im_x_tau8_re;
        re[i + 4 * ell] = tau1_re - tau3_re + tau5_re - tau7_re;
        im[i + 4 * ell] = tau1_im - tau3_im + tau5_im - tau7_im;
        re[i + 5 * ell] = tau2_re - a_re_x_tau4_re + a_im_x_tau4_im + tau6_im - b_re_x_tau8_re + b_im * tau8_im;
        im[i + 5 * ell] = tau2_im - a_re_x_tau4_im - a_im_x_tau4_re - tau6_re - b_re_x_tau8_im - b_im * tau8_re;
        re[i + 6 * ell] = tau1_re - tau3_im - tau5_re + tau7_im;
        im[i + 6 * ell] = tau1_im + tau3_re - tau5_im - tau7_re;
        re[i + 7 * ell] = tau2_re - b_re_x_tau4_re + b_im_x_tau4_im - tau6_im - a_re_x_tau8_re + a_im_x_tau8_im;
        im[i + 7 * ell] = tau2_im - b_re_x_tau4_im - b_im_x_tau4_re + tau6_re - a_re_x_tau8_im - a_im_x_tau8_re;
      }
    }
  }
}

void apply_conjugate_butterflies_r8(double* restrict re, double* restrict im, int n,
                                    const double** twiddle_re, const double** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 3; // Divide by 8
    t++;
  } // t = log8 n

  double tmp = 1. / sqrt(2.);
  double a_re =  tmp;
  double a_im = -tmp;
  double b_re = -tmp;
  double b_im = -tmp;
    
  for (int q = 0; q < t; ++q) {
    int k = 1 << (3 * (q + 1));  // k = 8^(q+1)
    int ell = k / 8;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double x1_re = re[i];
        double x1_im = im[i];
        double x2_re = re[i + ell];
        double x2_im = im[i + ell];
        double x3_re = re[i + 2 * ell];
        double x3_im = im[i + 2 * ell];
        double x4_re = re[i + 3 * ell];
        double x4_im = im[i + 3 * ell];
        double x5_re = re[i + 4 * ell];
        double x5_im = im[i + 4 * ell];
        double x6_re = re[i + 5 * ell];
        double x6_im = im[i + 5 * ell];
        double x7_re = re[i + 6 * ell];
        double x7_im = im[i + 6 * ell];
        double x8_re = re[i + 7 * ell];
        double x8_im = im[i + 7 * ell];

        double z1_re = x1_re;
        double z1_im = x1_im;
        double z2_re = twiddle_re[q][j] * x2_re + twiddle_im[q][j] * x2_im;
        double z2_im = twiddle_re[q][j] * x2_im - twiddle_im[q][j] * x2_re;
        double z3_re = twiddle_re[q][j + ell] * x3_re + twiddle_im[q][j + ell] * x3_im;
        double z3_im = twiddle_re[q][j + ell] * x3_im - twiddle_im[q][j + ell] * x3_re;
        double z4_re = twiddle_re[q][j + 2 * ell] * x4_re + twiddle_im[q][j + 2 * ell] * x4_im;
        double z4_im = twiddle_re[q][j + 2 * ell] * x4_im - twiddle_im[q][j + 2 * ell] * x4_re;
        double z5_re = twiddle_re[q][j + 3 * ell] * x5_re + twiddle_im[q][j + 3 * ell] * x5_im;
        double z5_im = twiddle_re[q][j + 3 * ell] * x5_im - twiddle_im[q][j + 3 * ell] * x5_re;
        double z6_re = twiddle_re[q][j + 4 * ell] * x6_re + twiddle_im[q][j + 4 * ell] * x6_im;
        double z6_im = twiddle_re[q][j + 4 * ell] * x6_im - twiddle_im[q][j + 4 * ell] * x6_re;
        double z7_re = twiddle_re[q][j + 5 * ell] * x7_re + twiddle_im[q][j + 5 * ell] * x7_im;
        double z7_im = twiddle_re[q][j + 5 * ell] * x7_im - twiddle_im[q][j + 5 * ell] * x7_re;
        double z8_re = twiddle_re[q][j + 6 * ell] * x8_re + twiddle_im[q][j + 6 * ell] * x8_im;
        double z8_im = twiddle_re[q][j + 6 * ell] * x8_im - twiddle_im[q][j + 6 * ell] * x8_re;

        double tau1_re = z1_re + z5_re;
        double tau1_im = z1_im + z5_im;
        double tau2_re = z1_re - z5_re;
        double tau2_im = z1_im - z5_im;
        double tau3_re = z2_re + z6_re;
        double tau3_im = z2_im + z6_im;
        double tau4_re = z2_re - z6_re;
        double tau4_im = z2_im - z6_im;
        double tau5_re = z3_re + z7_re;
        double tau5_im = z3_im + z7_im;
        double tau6_re = z3_re - z7_re;
        double tau6_im = z3_im - z7_im;
        double tau7_re = z4_re + z8_re;
        double tau7_im = z4_im + z8_im;
        double tau8_re = z4_re - z8_re;
        double tau8_im = z4_im - z8_im;

        re[i]           = tau1_re + tau3_re + tau5_re + tau7_re;
        im[i]           = tau1_im + tau3_im + tau5_im + tau7_im;        
        re[i + ell]     = tau2_re + a_re * tau4_re + a_im * tau4_im - tau6_im + b_re * tau8_re + b_im * tau8_im;
        im[i + ell]     = tau2_im + a_re * tau4_im - a_im * tau4_re + tau6_re + b_re * tau8_im - b_im * tau8_re;
        re[i + 2 * ell] = tau1_re - tau3_im - tau5_re + tau7_im;
        im[i + 2 * ell] = tau1_im + tau3_re - tau5_im - tau7_re;
        re[i + 3 * ell] = tau2_re + b_re * tau4_re + b_im * tau4_im + tau6_im + a_re * tau8_re + a_im * tau8_im;
        im[i + 3 * ell] = tau2_im + b_re * tau4_im - b_im * tau4_re - tau6_re + a_re * tau8_im - a_im * tau8_re;
        re[i + 4 * ell] = tau1_re - tau3_re + tau5_re - tau7_re;
        im[i + 4 * ell] = tau1_im - tau3_im + tau5_im - tau7_im;
        re[i + 5 * ell] = tau2_re - a_re * tau4_re - a_im * tau4_im - tau6_im - b_re * tau8_re - b_im * tau8_im;
        im[i + 5 * ell] = tau2_im - a_re * tau4_im + a_im * tau4_re + tau6_re - b_re * tau8_im + b_im * tau8_re;
        re[i + 6 * ell] = tau1_re + tau3_im - tau5_re - tau7_im;
        im[i + 6 * ell] = tau1_im - tau3_re - tau5_im + tau7_re;
        re[i + 7 * ell] = tau2_re - b_re * tau4_re - b_im * tau4_im + tau6_im - a_re * tau8_re - a_im * tau8_im;
        im[i + 7 * ell] = tau2_im - b_re * tau4_im + b_im * tau4_re - tau6_re - a_re * tau8_im + a_im * tau8_re;
      }
    }
  }
}

void apply_transposed_butterflies_r8(double* restrict re, double* restrict im, int n,
                                     const double** twiddle_re, const double** twiddle_im) {
  int t = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 3; // Divide by 8
    t++;
  } // t = log8 n

  double tmp = 1. / sqrt(2.);
  double a_re =  tmp;
  double a_im = -tmp;
  double b_re = -tmp;
  double b_im = -tmp;
    
  for (int q = 0; q < t; ++q) {
    int k = 1 << (3 * (q + 1));  // k = 8^(q+1)
    int ell = k / 8;
    int nb = n / k;
    for (int b = 0; b < nb; ++b) {
      int bk = b * k;

      for (int j = 0; j < ell; ++j) {
        int i = bk + j;

        double z1_re = re[i];
        double z1_im = im[i];
        double z2_re = re[i + ell];
        double z2_im = im[i + ell];
        double z3_re = re[i + 2 * ell];
        double z3_im = im[i + 2 * ell];
        double z4_re = re[i + 3 * ell];
        double z4_im = im[i + 3 * ell];
        double z5_re = re[i + 4 * ell];
        double z5_im = im[i + 4 * ell];
        double z6_re = re[i + 5 * ell];
        double z6_im = im[i + 5 * ell];
        double z7_re = re[i + 6 * ell];
        double z7_im = im[i + 6 * ell];
        double z8_re = re[i + 7 * ell];
        double z8_im = im[i + 7 * ell];

        double tau1_re = z1_re + z5_re;
        double tau1_im = z1_im + z5_im;
        double tau2_re = twiddle_re[q][j] * (z1_re - z5_re) - twiddle_im[q][j] * (z1_im - z5_im);
        double tau2_im = twiddle_re[q][j] * (z1_im - z5_im) + twiddle_im[q][j] * (z1_re - z5_re);
        double tau3_re = z2_re + z6_re;
        double tau3_im = z2_im + z6_im;
        double tau4_re = twiddle_re[q][j] * (z2_re - z6_re) - twiddle_im[q][j] * (z2_im - z6_im);
        double tau4_im = twiddle_re[q][j] * (z2_im - z6_im) + twiddle_im[q][j] * (z2_re - z6_re);
        double tau5_re = z3_re + z7_re;
        double tau5_im = z3_im + z7_im;
        double tau6_re = twiddle_re[q][j] * (z3_re - z7_re) - twiddle_im[q][j] * (z3_im - z7_im);
        double tau6_im = twiddle_re[q][j] * (z3_im - z7_im) + twiddle_im[q][j] * (z3_re - z7_re);
        double tau7_re = z4_re + z8_re;
        double tau7_im = z4_im + z8_im;
        double tau8_re = twiddle_re[q][j] * (z4_re - z8_re) - twiddle_im[q][j] * (z4_im - z8_im);
        double tau8_im = twiddle_re[q][j] * (z4_im - z8_im) + twiddle_im[q][j] * (z4_re - z8_re);

        re[i]           = tau1_re + tau3_re + tau5_re + tau7_re;
        re[i]           = tau1_im + tau3_im + tau5_im + tau7_im;
        re[i + ell]     = tau2_re + a_re * tau4_re - a_im * tau4_im + tau6_im + b_re * tau8_re - b_im * tau8_im;
        im[i + ell]     = tau2_im + a_re * tau4_im + a_im * tau4_re - tau6_re + b_re * tau8_im + b_im * tau8_re;
        re[i + 2 * ell] = twiddle_re[q][j + ell] * (tau1_re + tau3_im - tau5_re - tau7_im) - twiddle_im[q][j + ell] * (tau1_im - tau3_re - tau5_im + tau7_re);
        im[i + 2 * ell] = twiddle_re[q][j + ell] * (tau1_im - tau3_re - tau5_im + tau7_re) + twiddle_im[q][j + ell] * (tau1_re + tau3_im - tau5_re - tau7_im);
        re[i + 3 * ell] = twiddle_re[q][j + ell] * (tau2_re + b_re * tau4_re - b_im * tau4_im - tau6_im + a_re * tau8_re - a_im * tau8_im)
                        - twiddle_im[q][j + ell] * (tau2_im + b_re * tau4_im + b_im * tau4_re + tau6_re + a_re * tau8_im + a_im * tau8_re);
        im[i + 3 * ell] = twiddle_re[q][j + ell] * (tau2_im + b_re * tau4_im + b_im * tau4_re + tau6_re + a_re * tau8_im + a_im * tau8_re)
                        + twiddle_im[q][j + ell] * (tau2_re + b_re * tau4_re - b_im * tau4_im - tau6_im + a_re * tau8_re - a_im * tau8_im);
        re[i + 4 * ell] = twiddle_re[q][j + 3 * ell] * (tau1_re - tau3_re + tau5_re - tau7_re) - twiddle_im[q][j + 3 * ell] * (tau1_im - tau3_im + tau5_im - tau7_im);
        im[i + 4 * ell] = twiddle_re[q][j + 3 * ell] * (tau1_im - tau3_im + tau5_im - tau7_im) + twiddle_im[q][j + 3 * ell] * (tau1_re - tau3_re + tau5_re - tau7_re);
        re[i + 5 * ell] = twiddle_re[q][j + 3 * ell] * (tau2_re - a_re * tau4_re + a_im * tau4_im + tau6_im - b_re * tau8_re + b_im * tau8_im)
                        - twiddle_im[q][j + 3 * ell] * (tau2_im - a_re * tau4_im - a_im * tau4_re - tau6_re - b_re * tau8_im - b_im * tau8_re);
        im[i + 5 * ell] = twiddle_im[q][j + 3 * ell] * (tau2_re - a_re * tau4_re + a_im * tau4_im + tau6_im - b_re * tau8_re + b_im * tau8_im)
                        + twiddle_re[q][j + 3 * ell] * (tau2_im - a_re * tau4_im - a_im * tau4_re - tau6_re - b_re * tau8_im - b_im * tau8_re);
        re[i + 6 * ell] = twiddle_re[q][j + 5 * ell] * (tau1_re - tau3_im - tau5_re + tau7_im) - twiddle_im[q][j + 5 * ell] * (tau1_im + tau3_re - tau5_im - tau7_re);
        im[i + 6 * ell] = twiddle_im[q][j + 5 * ell] * (tau1_re - tau3_im - tau5_re + tau7_im) + twiddle_re[q][j + 5 * ell] * (tau1_im + tau3_re - tau5_im - tau7_re);
        re[i + 7 * ell] = twiddle_re[q][j + 5 * ell] * (tau2_re - b_re * tau4_re + b_im * tau4_im - tau6_im - a_re * tau8_re + a_im * tau8_im)
                        - twiddle_im[q][j + 5 * ell] * (tau2_im - b_re * tau4_im - b_im * tau4_re + tau6_re - a_re * tau8_im - a_im * tau8_re);
        im[i + 7 * ell] = twiddle_im[q][j + 5 * ell] * (tau2_re - b_re * tau4_re + b_im * tau4_im - tau6_im - a_re * tau8_re + a_im * tau8_im)
                        + twiddle_re[q][j + 5 * ell] * (tau2_im - b_re * tau4_im - b_im * tau4_re + tau6_re - a_re * tau8_im - a_im * tau8_re);
      }
    }
  }
}