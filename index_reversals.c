#include <stdlib.h>
#include <complex.h>

int reverse_index_r2(int x, int t) {
    int r = 0;
    for (int i = 0; i < t; ++i) {
        r = (r << 1) | (x & 1);  // Extract lowest 1 bit and shift result
        x >>= 1;  // Remove lowest 1 bit from x
    }
    return r;
}

int* precompute_index_reversal_permutation_r2(int n) {
    int t = __builtin_ctz(n); // t = log2 n
    int* rho = (int*) malloc(n * sizeof(int));
    for (int i = 0; i < n; ++i)
        rho[i] = reverse_index_r2(i, t);
    return rho;
}

int reverse_index_r4(int x, int t) {
    int r = 0;
    for (int i = 0; i < t; ++i) {
        r = (r << 2) | (x & 3);  // Extract lowest 2 bits and shift result
        x >>= 2;  // Remove lowest 2 bits from x
    }
    return r;
}

int* precompute_index_reversal_permutation_r4(int n) {
    int t = 0;
    int temp = n;
    while (temp > 1) {
        temp >>= 2;  // Divide by 4
        t++;
    } // t = log4 n
    int* rho = (int*) malloc(n * sizeof(int));
    for (int i = 0; i < n; ++i)
        rho[i] = reverse_index_r4(i, t);
    return rho;
}

int reverse_index_r8(int x, int t) {
    int r = 0;
    for (int i = 0; i < t; ++i) {
        r = (r << 3) | (x & 7);  // Extract lowest 3 bits and shift result
        x >>= 3;  // Remove lowest 3 bits from x
    }
    return r;
}

int* precompute_index_reversal_permutation_r8(int n) {
    int t = 0;
    int temp = n;
    while (temp > 1) {
        temp >>= 3;  // Divide by 8
        t++;
    } // t = log8 n
    int* rho = (int*) malloc(n * sizeof(int));
    for (int i = 0; i < n; ++i)
        rho[i] = reverse_index_r8(i, t);
    return rho;
}

void apply_permutation(double complex* x, int n, const int* rho) {
    for (int j = 0; j < n; ++j) {
        int r = rho[j];
        if (j < r) {
            double complex tmp = x[j];
            x[j] = x[r];
            x[r] = tmp;
        }
    }
}