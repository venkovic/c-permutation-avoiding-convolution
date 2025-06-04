# Copyright (c) 2025 Nicolas Venkovic
#
# This file is part of c-permutation-avoiding-convolution.
# 
# This file is licensed under the MIT License.
# For the full license text, see the LICENSE file in the root directory.

# Compiler and flags
CC = gcc
CFLAGS = -O3 -ffast-math -funroll-loops -Wall -Wextra
LDFLAGS = -lfftw3 -lm

# Common source files (without the main files)
COMMON_SRC = twiddle_factors.c butterflies.c butterflies_2d.c butterflies_3d.c index_reversals.c index_reversals_2d.c index_reversals_3d.c fft.c fft_2d.c fft_3d.c convolution.c convolution_2d.c convolution_3d.c misc.c
COMMON_OBJ = $(COMMON_SRC:.c=.o)
DEPS = fft.h

# Executables
TARGET1 = bench1d_forward_fft
TARGET2 = bench1d_backward_fft
TARGET3 = bench1d_convolution
TARGET4 = bench2d_forward_fft
TARGET5 = bench2d_backward_fft
TARGET6 = bench2d_convolution
TARGET7 = bench3d_forward_fft
TARGET8 = bench3d_backward_fft
TARGET9 = bench3d_convolution
TARGETS = $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8) $(TARGET9)

# Default rule - build both executables
all: $(TARGETS)

# Build forward FFT benchmark
$(TARGET1): $(COMMON_OBJ) bench1d_forward_fft.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build backward FFT benchmark
$(TARGET2): $(COMMON_OBJ) bench1d_backward_fft.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build convolution benchmark
$(TARGET3): $(COMMON_OBJ) bench1d_convolution.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build 2D forward FFT benchmark
$(TARGET4): $(COMMON_OBJ) bench2d_forward_fft.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build 2D backward FFT benchmark
$(TARGET5): $(COMMON_OBJ) bench2d_backward_fft.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build 2D convolution benchmark
$(TARGET6): $(COMMON_OBJ) bench2d_convolution.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build 3D forward FFT benchmark
$(TARGET7): $(COMMON_OBJ) bench3d_forward_fft.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build 3D backward FFT benchmark
$(TARGET8): $(COMMON_OBJ) bench3d_backward_fft.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build 3D convolution benchmark
$(TARGET9): $(COMMON_OBJ) bench3d_convolution.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile .c to .o
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(COMMON_OBJ) bench1d_forward_fft.o bench1d_backward_fft.o bench1d_convolution.o bench2d_forward_fft.o bench2d_backward_fft.o bench2d_convolution.o bench3d_forward_fft.o bench3d_backward_fft.o bench3d_convolution.o $(TARGETS)

# Run benchmarks
run-forward: $(TARGET1)
	./$(TARGET1)

run-backward: $(TARGET2)
	./$(TARGET2)

run-convolution: $(TARGET3)
	./$(TARGET3)

run-forward-2d: $(TARGET4)
	./$(TARGET4)

run-backward-2d: $(TARGET5)
	./$(TARGET5)

run-convolution-2d: $(TARGET6)
	./$(TARGET6)

run-forward-3d: $(TARGET7)
	./$(TARGET7)

run-backward-3d: $(TARGET8)
	./$(TARGET8)

run-convolution-3d: $(TARGET9)
	./$(TARGET9)

run-all: $(TARGETS)
	./$(TARGET1)
	./$(TARGET2)
	./$(TARGET3)
	./$(TARGET4)
	./$(TARGET5)
	./$(TARGET6)
	./$(TARGET7)
	./$(TARGET8)
	./$(TARGET9)

# Individual targets (useful for building just one)
forward: $(TARGET1)
backward: $(TARGET2)
convolution: $(TARGET3)
forward-2d: $(TARGET4)
backward-2d: $(TARGET5)
convolution-2d: $(TARGET6)
forward-3d: $(TARGET7)
backward-3d: $(TARGET8)
convolution-3d: $(TARGET9)

.PHONY: all clean run-forward run-backward run-convolution run-forward-2d run-backward-2d run-convolution-2d run-forward-3d run-backward-3d run-convolution-3d run-all forward backward convolution forward-2d backward-2d convolution-2d forward-3d backward-3d convolution-3d