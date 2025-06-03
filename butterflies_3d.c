#include "fft.h"

void apply_butterflies_r2_3d(double* restrict re, double* restrict im, 
                             int n1, int n2, int n3,
                             const double** twiddle_re1, const double** twiddle_im1,
                             const double** twiddle_re2, const double** twiddle_im2,
                             const double** twiddle_re3, const double** twiddle_im3) {
  
  // Temporary arrays for processing
  double* temp_re = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  double* temp_im = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  
  // Step 1: Apply 1D FFT along dimension 3 (innermost, contiguous)
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      double* re_slice = re + (i1 * n2 + i2) * n3;
      double* im_slice = im + (i1 * n2 + i2) * n3;
      
      apply_butterflies_r2(re_slice, im_slice, n3, twiddle_re1, twiddle_im1);
    }
  }
  
  // Step 2: Apply 1D FFT along dimension 2
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i3 = 0; i3 < n3; ++i3) {
      // Extract data along dimension 2
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i2] = re[idx];
        temp_im[i2] = im[idx];
      }
      
      // Apply 1D FFT
      apply_butterflies_r2(temp_re, temp_im, n2, twiddle_re2, twiddle_im2);
      
      // Copy back
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i2];
        im[idx] = temp_im[i2];
      }
    }
  }
  
  // Step 3: Apply 1D FFT along dimension 1 (outermost)
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i3 = 0; i3 < n3; ++i3) {
      // Extract data along dimension 1
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i1] = re[idx];
        temp_im[i1] = im[idx];
      }
      
      // Apply 1D FFT
      apply_butterflies_r2(temp_re, temp_im, n1, twiddle_re3, twiddle_im3);
      
      // Copy back
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i1];
        im[idx] = temp_im[i1];
      }
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_conjugate_butterflies_r2_3d(double* restrict re, double* restrict im, 
                                       int n1, int n2, int n3,
                                       const double** twiddle_re1, const double** twiddle_im1,
                                       const double** twiddle_re2, const double** twiddle_im2,
                                       const double** twiddle_re3, const double** twiddle_im3) {
  
  double* temp_re = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  double* temp_im = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  
  // Step 1: Apply 1D conjugate FFT along dimension 3
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      double* re_slice = re + (i1 * n2 + i2) * n3;
      double* im_slice = im + (i1 * n2 + i2) * n3;
      
      apply_conjugate_butterflies_r2(re_slice, im_slice, n3, twiddle_re1, twiddle_im1);
    }
  }
  
  // Step 2: Apply 1D conjugate FFT along dimension 2
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i2] = re[idx];
        temp_im[i2] = im[idx];
      }
      
      apply_conjugate_butterflies_r2(temp_re, temp_im, n2, twiddle_re2, twiddle_im2);
      
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i2];
        im[idx] = temp_im[i2];
      }
    }
  }
  
  // Step 3: Apply 1D conjugate FFT along dimension 1
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i1] = re[idx];
        temp_im[i1] = im[idx];
      }
      
      apply_conjugate_butterflies_r2(temp_re, temp_im, n1, twiddle_re3, twiddle_im3);
      
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i1];
        im[idx] = temp_im[i1];
      }
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_transposed_butterflies_r2_3d(double* restrict re, double* restrict im, 
                                        int n1, int n2, int n3,
                                        const double** twiddle_re1, const double** twiddle_im1,
                                        const double** twiddle_re2, const double** twiddle_im2,
                                        const double** twiddle_re3, const double** twiddle_im3) {
  
  double* temp_re = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  double* temp_im = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  
  // Step 1: Apply 1D transposed FFT along dimension 3
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      double* re_slice = re + (i1 * n2 + i2) * n3;
      double* im_slice = im + (i1 * n2 + i2) * n3;
      
      apply_transposed_butterflies_r2(re_slice, im_slice, n3, twiddle_re1, twiddle_im1);
    }
  }
  
  // Step 2: Apply 1D transposed FFT along dimension 2
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i2] = re[idx];
        temp_im[i2] = im[idx];
      }
      
      apply_transposed_butterflies_r2(temp_re, temp_im, n2, twiddle_re2, twiddle_im2);
      
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i2];
        im[idx] = temp_im[i2];
      }
    }
  }
  
  // Step 3: Apply 1D transposed FFT along dimension 1
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i1] = re[idx];
        temp_im[i1] = im[idx];
      }
      
      apply_transposed_butterflies_r2(temp_re, temp_im, n1, twiddle_re3, twiddle_im3);
      
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i1];
        im[idx] = temp_im[i1];
      }
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_butterflies_r4_3d(double* restrict re, double* restrict im, 
                             int n1, int n2, int n3,
                             const double** twiddle_re1, const double** twiddle_im1,
                             const double** twiddle_re2, const double** twiddle_im2,
                             const double** twiddle_re3, const double** twiddle_im3) {
  
  double* temp_re = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  double* temp_im = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  
  // Step 1: Apply 1D radix-4 FFT along dimension 3
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      double* re_slice = re + (i1 * n2 + i2) * n3;
      double* im_slice = im + (i1 * n2 + i2) * n3;
      
      apply_butterflies_r4(re_slice, im_slice, n3, twiddle_re1, twiddle_im1);
    }
  }
  
  // Step 2: Apply 1D radix-4 FFT along dimension 2
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i2] = re[idx];
        temp_im[i2] = im[idx];
      }
      
      apply_butterflies_r4(temp_re, temp_im, n2, twiddle_re2, twiddle_im2);
      
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i2];
        im[idx] = temp_im[i2];
      }
    }
  }
  
  // Step 3: Apply 1D radix-4 FFT along dimension 1
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i1] = re[idx];
        temp_im[i1] = im[idx];
      }
      
      apply_butterflies_r4(temp_re, temp_im, n1, twiddle_re3, twiddle_im3);
      
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i1];
        im[idx] = temp_im[i1];
      }
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_conjugate_butterflies_r4_3d(double* restrict re, double* restrict im, 
                                       int n1, int n2, int n3,
                                       const double** twiddle_re1, const double** twiddle_im1,
                                       const double** twiddle_re2, const double** twiddle_im2,
                                       const double** twiddle_re3, const double** twiddle_im3) {
  
  double* temp_re = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  double* temp_im = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  
  // Step 1: Apply 1D conjugate radix-4 FFT along dimension 3
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      double* re_slice = re + (i1 * n2 + i2) * n3;
      double* im_slice = im + (i1 * n2 + i2) * n3;
      
      apply_conjugate_butterflies_r4(re_slice, im_slice, n3, twiddle_re1, twiddle_im1);
    }
  }
  
  // Step 2: Apply 1D conjugate radix-4 FFT along dimension 2
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i2] = re[idx];
        temp_im[i2] = im[idx];
      }
      
      apply_conjugate_butterflies_r4(temp_re, temp_im, n2, twiddle_re2, twiddle_im2);
      
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i2];
        im[idx] = temp_im[i2];
      }
    }
  }
  
  // Step 3: Apply 1D conjugate radix-4 FFT along dimension 1
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i1] = re[idx];
        temp_im[i1] = im[idx];
      }
      
      apply_conjugate_butterflies_r4(temp_re, temp_im, n1, twiddle_re3, twiddle_im3);
      
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i1];
        im[idx] = temp_im[i1];
      }
    }
  }
  
  free(temp_re);
  free(temp_im);
}

void apply_transposed_butterflies_r4_3d(double* restrict re, double* restrict im, 
                                        int n1, int n2, int n3,
                                        const double** twiddle_re1, const double** twiddle_im1,
                                        const double** twiddle_re2, const double** twiddle_im2,
                                        const double** twiddle_re3, const double** twiddle_im3) {
  
  double* temp_re = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  double* temp_im = (double*)malloc((n1 > n2 ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3)) * sizeof(double));
  
  // Step 1: Apply 1D transposed radix-4 FFT along dimension 3
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      double* re_slice = re + (i1 * n2 + i2) * n3;
      double* im_slice = im + (i1 * n2 + i2) * n3;
      
      apply_transposed_butterflies_r4(re_slice, im_slice, n3, twiddle_re1, twiddle_im1);
    }
  }
  
  // Step 2: Apply 1D transposed radix-4 FFT along dimension 2
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i2] = re[idx];
        temp_im[i2] = im[idx];
      }
      
      apply_transposed_butterflies_r4(temp_re, temp_im, n2, twiddle_re2, twiddle_im2);
      
      for (int i2 = 0; i2 < n2; ++i2) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i2];
        im[idx] = temp_im[i2];
      }
    }
  }
  
  // Step 3: Apply 1D transposed radix-4 FFT along dimension 1
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        temp_re[i1] = re[idx];
        temp_im[i1] = im[idx];
      }
      
      apply_transposed_butterflies_r4(temp_re, temp_im, n1, twiddle_re3, twiddle_im3);
      
      for (int i1 = 0; i1 < n1; ++i1) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        re[idx] = temp_re[i1];
        im[idx] = temp_im[i1];
      }
    }
  }
  
  free(temp_re);
  free(temp_im);
}