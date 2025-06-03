#include "fft.h"

void apply_butterflies_r2_2d(double* restrict re, double* restrict im, 
                             int n1, int n2,
                             const double** twiddle_re1, const double** twiddle_im1,
                             const double** twiddle_re2, const double** twiddle_im2) {
  
  // Temporary arrays for column processing
  double* temp_re = (double*)malloc(n1 * sizeof(double));
  double* temp_im = (double*)malloc(n2 * sizeof(double));
  
  // Step 1: Apply 1D butterflies to each row
  for (int i1 = 0; i1 < n1; ++i1) {
    double* re1 = re + i1 * n2;  // Point to start of current row
    double* im1 = im + i1 * n2;
    
    apply_butterflies_r2(re1, im1, n2, twiddle_re1, twiddle_im1);
  }
  
  // Step 2: Apply 1D butterflies to each column
  for (int i2 = 0; i2 < n2; ++i2) {
    // Extract column data into temporary arrays
    for (int i1 = 0; i1 < n1; ++i1) {
      temp_re[i1] = re[i1 * n2 + i2];
      temp_im[i1] = im[i1 * n2 + i2];
    }
    
    // Apply 1D butterflies to the column
    apply_butterflies_r2(temp_re, temp_im, n1, twiddle_re2, twiddle_im2);
    
    // Copy results back to the 2D array
    for (int i1 = 0; i1 < n1; ++i1) {
      re[i1 * n2 + i2] = temp_re[i1];
      im[i1 * n2 + i2] = temp_im[i1];
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_conjugate_butterflies_r2_2d(double* restrict re, double* restrict im, 
                                       int n1, int n2,
                                       const double** twiddle_re1, const double** twiddle_im1,
                                       const double** twiddle_re2, const double** twiddle_im2) {
  
  // Temporary arrays for column processing
  double* temp_re = (double*)malloc(n1 * sizeof(double));
  double* temp_im = (double*)malloc(n1 * sizeof(double));
  
  // Step 1: Apply 1D conjugate butterflies to each row
  for (int i1 = 0; i1 < n1; ++i1) {
    double* re1 = re + i1 * n2;  // Point to start of current row
    double* im1 = im + i1 * n2;
    
    apply_conjugate_butterflies_r2(re1, im1, n2, twiddle_re1, twiddle_im1);
  }
  
  // Step 2: Apply 1D conjugate butterflies to each column
  for (int i2 = 0; i2 < n2; ++i2) {
    // Extract column data into temporary arrays
    for (int i1 = 0; i1 < n1; ++i1) {
      temp_re[i1] = re[i1 * n2 + i2];
      temp_im[i1] = im[i1 * n2 + i2];
    }
    
    // Apply 1D conjugate butterflies to the column
    apply_conjugate_butterflies_r2(temp_re, temp_im, n1, twiddle_re2, twiddle_im2);
    
    // Copy results back to the 2D array
    for (int i1 = 0; i1 < n1; ++i1) {
      re[i1 * n2 + i2] = temp_re[i1];
      im[i1 * n2 + i2] = temp_im[i1];
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_transposed_butterflies_r2_2d(double* restrict re, double* restrict im, 
                                        int n1, int n2,
                                        const double** twiddle_re1, const double** twiddle_im1,
                                        const double** twiddle_re2, const double** twiddle_im2) {
  
  // Temporary arrays for column processing
  double* temp_re = (double*)malloc(n1 * sizeof(double));
  double* temp_im = (double*)malloc(n1 * sizeof(double));
  
  // Step 1: Apply 1D transposed butterflies to each row
  for (int i1 = 0; i1 < n1; ++i1) {
    double* re1 = re + i1 * n2;  // Point to start of current row
    double* im1 = im + i1 * n2;
    
    apply_transposed_butterflies_r2(re1, im1, n2, twiddle_re1, twiddle_im1);
  }
  
  // Step 2: Apply 1D transposed butterflies to each column
  for (int i2 = 0; i2 < n2; ++i2) {
    // Extract column data into temporary arrays
    for (int i1 = 0; i1 < n1; ++i1) {
      temp_re[i1] = re[i1 * n2 + i2];
      temp_im[i1] = im[i1 * n2 + i2];
    }
    
    // Apply 1D transposed butterflies to the column
    apply_transposed_butterflies_r2(temp_re, temp_im, n1, twiddle_re2, twiddle_im2);
    
    // Copy results back to the 2D array
    for (int i1 = 0; i1 < n1; ++i1) {
      re[i1 * n2 + i2] = temp_re[i1];
      im[i1 * n2 + i2] = temp_im[i1];
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_butterflies_r4_2d(double* restrict re, double* restrict im, 
                             int n1, int n2,
                             const double** twiddle_re1, const double** twiddle_im1,
                             const double** twiddle_re2, const double** twiddle_im2) {
  
  // Temporary arrays for column processing
  double* temp_re = (double*)malloc(n1 * sizeof(double));
  double* temp_im = (double*)malloc(n1 * sizeof(double));
  
  // Step 1: Apply 1D radix-4 butterflies to each row
  for (int i1 = 0; i1 < n1; ++i1) {
    double* re1 = re + i1 * n2;  // Point to start of current row
    double* im1 = im + i1 * n2;
    
    apply_butterflies_r4(re1, im1, n2, twiddle_re1, twiddle_im1);
  }
  
  // Step 2: Apply 1D radix-4 butterflies to each column
  for (int i2 = 0; i2 < n2; ++i2) {
    // Extract column data into temporary arrays
    for (int i1 = 0; i1 < n1; ++i1) {
      temp_re[i1] = re[i1 * n2 + i2];
      temp_im[i1] = im[i1 * n2 + i2];
    }
    
    // Apply 1D radix-4 butterflies to the column
    apply_butterflies_r4(temp_re, temp_im, n1, twiddle_re2, twiddle_im2);
    
    // Copy results back to the 2D array
    for (int i1 = 0; i1 < n1; ++i1) {
      re[i1 * n2 + i2] = temp_re[i1];
      im[i1 * n2 + i2] = temp_im[i1];
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_conjugate_butterflies_r4_2d(double* restrict re, double* restrict im, 
                                       int n1, int n2,
                                       const double** twiddle_re1, const double** twiddle_im1,
                                       const double** twiddle_re2, const double** twiddle_im2) {
  
  // Temporary arrays for column processing
  double* temp_re = (double*)malloc(n1 * sizeof(double));
  double* temp_im = (double*)malloc(n1 * sizeof(double));
  
  // Step 1: Apply 1D conjugate radix-4 butterflies to each row
  for (int i1 = 0; i1 < n1; ++i1) {
    double* re1 = re + i1 * n2;  // Point to start of current row
    double* im1 = im + i1 * n2;
    
    apply_conjugate_butterflies_r4(re1, im1, n2, twiddle_re1, twiddle_im1);
  }
  
  // Step 2: Apply 1D conjugate radix-4 butterflies to each column
  for (int i2 = 0; i2 < n2; ++i2) {
    // Extract column data into temporary arrays
    for (int i1 = 0; i1 < n1; ++i1) {
      temp_re[i1] = re[i1 * n2 + i2];
      temp_im[i1] = im[i1 * n2 + i2];
    }
    
    // Apply 1D conjugate radix-4 butterflies to the column
    apply_conjugate_butterflies_r4(temp_re, temp_im, n1, twiddle_re2, twiddle_im2);
    
    // Copy results back to the 2D array
    for (int i1 = 0; i1 < n1; ++i1) {
      re[i1 * n2 + i2] = temp_re[i1];
      im[i1 * n2 + i2] = temp_im[i1];
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_transposed_butterflies_r4_2d(double* restrict re, double* restrict im, 
                                        int n1, int n2,
                                        const double** twiddle_re1, const double** twiddle_im1,
                                        const double** twiddle_re2, const double** twiddle_im2) {
  
  // Temporary arrays for column processing
  double* temp_re = (double*)malloc(n1 * sizeof(double));
  double* temp_im = (double*)malloc(n1 * sizeof(double));
  
  // Step 1: Apply 1D transposed radix-4 butterflies to each row
  for (int i1 = 0; i1 < n1; ++i1) {
    double* re1 = re + i1 * n2;  // Point to start of current row
    double* im1 = im + i1 * n2;
    
    apply_transposed_butterflies_r4(re1, im1, n2, twiddle_re1, twiddle_im1);
  }
  
  // Step 2: Apply 1D transposed radix-4 butterflies to each column
  for (int i2 = 0; i2 < n2; ++i2) {
    // Extract column data into temporary arrays
    for (int i1 = 0; i1 < n1; ++i1) {
      temp_re[i1] = re[i1 * n2 + i2];
      temp_im[i1] = im[i1 * n2 + i2];
    }
    
    // Apply 1D transposed radix-4 butterflies FFT to the column
    apply_transposed_butterflies_r4(temp_re, temp_im, n1, twiddle_re2, twiddle_im2);
    
    // Copy results back to the 2D array
    for (int i1 = 0; i1 < n1; ++i1) {
      re[i1 * n2 + i2] = temp_re[i1];
      im[i1 * n2 + i2] = temp_im[i1];
    }
  }
  
  free(temp_re);
  free(temp_im);
}