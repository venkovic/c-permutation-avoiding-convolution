#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define ARRAY_SIZE (256 * 1024 * 1024)  // 256 MB - larger than L3 cache
#define ITERATIONS 5
#define WARMUP_ITERATIONS 2

double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec / 1e9;
}

// Sequential read test
double test_read_bandwidth(double *array, size_t size) {
  double start, end;
  double sum = 0.0;
    
  start = get_time();
  for (size_t i = 0; i < size; i++) {
    sum += array[i];
  }
  end = get_time();
    
  // Prevent compiler optimization
  if (sum == 0.0) printf(" ");
    
  return (size * sizeof(double)) / (end - start) / (1024.0 * 1024.0 * 1024.0);
}

// Sequential write test
double test_write_bandwidth(double *array, size_t size) {
  double start, end;
    
  start = get_time();
  for (size_t i = 0; i < size; i++) {
    array[i] = 1.0;
  }
  end = get_time();
    
  return (size * sizeof(double)) / (end - start) / (1024.0 * 1024.0 * 1024.0);
}

// Copy test (read + write)
double test_copy_bandwidth(double *src, double *dst, size_t size) {
  double start, end;
    
  start = get_time();
  for (size_t i = 0; i < size; i++) {
    dst[i] = src[i];
  }
  end = get_time();
    
  return (2 * size * sizeof(double)) / (end - start) / (1024.0 * 1024.0 * 1024.0);
}

int main() {
  printf("Single Core Memory Bandwidth Test\n");
  printf("=================================\n");
    
  // Get system info
  printf("Array size: %d MB\n", (int)(ARRAY_SIZE * sizeof(double) / (1024 * 1024)));
  printf("Iterations: %d (plus %d warmup)\n\n", ITERATIONS, WARMUP_ITERATIONS);
    
  // Allocate memory
  double *array1 = malloc(ARRAY_SIZE * sizeof(double));
  double *array2 = malloc(ARRAY_SIZE * sizeof(double));
    
  if (!array1 || !array2) {
    printf("Memory allocation failed!\n");
    return 1;
  }
    
  // Initialize arrays
  for (size_t i = 0; i < ARRAY_SIZE; i++) {
    array1[i] = (double)i;
    array2[i] = 0.0;
  }
    
  printf("Running tests...\n");
    
  // Warmup
  for (int i = 0; i < WARMUP_ITERATIONS; i++) {
    test_read_bandwidth(array1, ARRAY_SIZE);
    test_write_bandwidth(array1, ARRAY_SIZE);
    test_copy_bandwidth(array1, array2, ARRAY_SIZE);
  }
    
  // Read bandwidth test
  double read_total = 0.0;
  printf("\nRead Bandwidth:\n");
  for (int i = 0; i < ITERATIONS; i++) {
    double bw = test_read_bandwidth(array1, ARRAY_SIZE);
    printf("  Run %d: %.2f GB/s\n", i+1, bw);
    read_total += bw;
  }
  double read_avg = read_total / ITERATIONS;
    
  // Write bandwidth test
  double write_total = 0.0;
  printf("\nWrite Bandwidth:\n");
  for (int i = 0; i < ITERATIONS; i++) {
    double bw = test_write_bandwidth(array1, ARRAY_SIZE);
    printf("  Run %d: %.2f GB/s\n", i+1, bw);
    write_total += bw;
  }
  double write_avg = write_total / ITERATIONS;
    
  // Copy bandwidth test
  double copy_total = 0.0;
  printf("\nCopy Bandwidth:\n");
  for (int i = 0; i < ITERATIONS; i++) {
    double bw = test_copy_bandwidth(array1, array2, ARRAY_SIZE);
    printf("  Run %d: %.2f GB/s\n", i+1, bw);
    copy_total += bw;
  }
  double copy_avg = copy_total / ITERATIONS;
    
  // Results summary
  printf("\n=================================\n");
  printf("RESULTS SUMMARY:\n");
  printf("=================================\n");
  printf("Average Read Bandwidth:  %.2f GB/s\n", read_avg);
  printf("Average Write Bandwidth: %.2f GB/s\n", write_avg);
  printf("Average Copy Bandwidth:  %.2f GB/s\n", copy_avg);
  printf("Peak Single-Core BW:    %.2f GB/s\n", 
  (read_avg > write_avg) ? read_avg : write_avg);
    
  free(array1);
  free(array2);
    
  return 0;
}