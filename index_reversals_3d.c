/*
 * Copyright (c) 2025 Nicolas Venkovic
 * 
 * This file is part of c-permutation-avoiding-convolution.
 * 
 * This file is licensed under the MIT License.
 * For the full license text, see the LICENSE file in the root directory of this project.
 */

#include <stdlib.h>

#include "fft.h"

int* precompute_index_reversal_permutation_r2_3d(int n1, int n2, int n3) {
  int n = n1 * n2 * n3;
  int* rho = (int*)malloc(n * sizeof(int));
  
  int t1 = __builtin_ctz(n1);
  int t2 = __builtin_ctz(n2);
  int t3 = __builtin_ctz(n3);
  
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      for (int i3 = 0; i3 < n3; ++i3) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        
        // Apply bit-reversal to all three dimensions independently
        int rev_i1 = reverse_index_r2(i1, t1);
        int rev_i2 = reverse_index_r2(i2, t2);
        int rev_i3 = reverse_index_r2(i3, t3);
        int rev_idx = (rev_i1 * n2 + rev_i2) * n3 + rev_i3;
        
        rho[idx] = rev_idx;
      }
    }
  }
  
  return rho;
}

int* precompute_index_reversal_permutation_r4_3d(int n1, int n2, int n3) {
  int n = n1 * n2 * n3;
  int* rho = (int*)malloc(n * sizeof(int));
  
  // Calculate log4 for all three dimensions
  int t1 = 0, temp1 = n1;
  while (temp1 > 1) { temp1 >>= 2; t1++; }  // t1 = log4(n1)
  
  int t2 = 0, temp2 = n2;
  while (temp2 > 1) { temp2 >>= 2; t2++; }  // t2 = log4(n2)
  
  int t3 = 0, temp3 = n3;
  while (temp3 > 1) { temp3 >>= 2; t3++; }  // t3 = log4(n3)
  
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      for (int i3 = 0; i3 < n3; ++i3) {
        int idx = (i1 * n2 + i2) * n3 + i3;
        
        // Apply radix-4 bit-reversal to all three dimensions independently
        int rev_i1 = reverse_index_r4(i1, t1);
        int rev_i2 = reverse_index_r4(i2, t2);
        int rev_i3 = reverse_index_r4(i3, t3);
        int rev_idx = (rev_i1 * n2 + rev_i2) * n3 + rev_i3;
        
        rho[idx] = rev_idx;
      }
    }
  }
  
  return rho;
}