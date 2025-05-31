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

                double x0_re = re[i];
                double x0_im = im[i];
                double x1_re = re[i + ell];
                double x1_im = im[i + ell];
                double x2_re = re[i + 2 * ell];
                double x2_im = im[i + 2 * ell];
                double x3_re = re[i + 3 * ell];
                double x3_im = im[i + 3 * ell];

                double alpha1_re = twiddle_re[q][j] * x1_re - twiddle_im[q][j] * x1_im;
                double alpha1_im = twiddle_re[q][j] * x1_im + twiddle_im[q][j] * x1_re;

                double alpha2_re = twiddle_re[q][j + ell] * x2_re - twiddle_im[q][j + ell] * x2_im;
                double alpha2_im = twiddle_re[q][j + ell] * x2_im + twiddle_im[q][j + ell] * x2_re;

                double alpha3_re = twiddle_re[q][j + 2 * ell] * x3_re - twiddle_im[q][j + 2 * ell] * x3_im;
                double alpha3_im = twiddle_re[q][j + 2 * ell] * x3_im + twiddle_im[q][j + 2 * ell] * x3_re;

                double tau0_re = x0_re + alpha2_re;
                double tau0_im = x0_im + alpha2_im;
                double tau1_re = x0_re - alpha2_re;
                double tau1_im = x0_im - alpha2_im;
                double tau2_re = alpha1_re + alpha3_re;
                double tau2_im = alpha1_im + alpha3_im;
                double tau3_re = alpha1_re - alpha3_re;
                double tau3_im = alpha1_im - alpha3_im;

                re[i]           = tau0_re + tau2_re;
                im[i]           = tau0_im + tau2_im;
                re[i + ell]     = tau1_re + tau3_im;
                im[i + ell]     = tau1_im - tau3_re;
                re[i + 2 * ell] = tau0_re - tau2_re;
                im[i + 2 * ell] = tau0_im - tau2_im;
                re[i + 3 * ell] = tau1_re - tau3_im;
                im[i + 3 * ell] = tau1_im + tau3_re;
            }
        }
    }
}

