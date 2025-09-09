# Copyright (c) 2025 Nicolas Venkovic
#
# This file is part of c-permutation-avoiding-convolution.
# 
# This file is licensed under the MIT License.
# For the full license text, see the LICENSE file in the root directory.

# Compiler and flags
CC = gcc
CFLAGS = -O3 -ffast-math -funroll-loops -Wall -Wextra

# FFTW configuration - try multiple sources
FFTW_CFLAGS := $(shell pkg-config --cflags fftw3 2>/dev/null || \
                  echo "-I$(shell spack location -i fftw 2>/dev/null)/include" || \
                  echo "-I./local/include")
FFTW_LIBS := $(shell pkg-config --libs fftw3 2>/dev/null || \
               echo "-L$(shell spack location -i fftw 2>/dev/null)/lib -lfftw3" || \
               echo "-L./local/lib -lfftw3")

# Add FFTW flags to compiler flags
CFLAGS += $(FFTW_CFLAGS)
LDFLAGS = $(FFTW_LIBS) -lm

# Alternative: Static linking (uncomment to use)
# LDFLAGS = -L$(shell spack location -i fftw)/lib -Wl,-Bstatic -lfftw3 -Wl,-Bdynamic -lm

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
TARGET10 = bench2d_convolution_anisotropic
TARGET11 = bench3d_convolution_anisotropic
MEMORY_BW_TEST = memory_bw_test
TARGETS = $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8) $(TARGET9) $(TARGET10) $(TARGET11) $(MEMORY_BW_TEST)

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

# Build 2D anisotropic convolution benchmark
$(TARGET10): $(COMMON_OBJ) bench2d_convolution_anisotropic.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build 3D anisotropic convolution benchmark
$(TARGET11): $(COMMON_OBJ) bench3d_convolution_anisotropic.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build memory bandwidth test (no FFTW dependency)
$(MEMORY_BW_TEST): memory_bw_test.o
	$(CC) -O3 -o $@ $< -lm

# Compile .c to .o
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

# Special rule for memory bandwidth test (no FFTW headers needed)
memory_bw_test.o: memory_bw_test.c
	$(CC) -O3 -Wall -Wextra -c $< -o $@

# Clean rule
clean:
	rm -f $(COMMON_OBJ) bench1d_forward_fft.o bench1d_backward_fft.o bench1d_convolution.o bench2d_forward_fft.o bench2d_backward_fft.o bench2d_convolution.o bench3d_forward_fft.o bench3d_backward_fft.o bench3d_convolution.o bench2d_convolution_anisotropic.o bench3d_convolution_anisotropic.o memory_bw_test.o $(TARGETS)

# Debug rule to check FFTW configuration
debug-fftw:
	@echo "FFTW_CFLAGS: $(FFTW_CFLAGS)"
	@echo "FFTW_LIBS: $(FFTW_LIBS)"
	@echo "Final CFLAGS: $(CFLAGS)"
	@echo "Final LDFLAGS: $(LDFLAGS)"

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

run-convolution-anisotropic-2d: $(TARGET10)
	./$(TARGET10)

run-convolution-anisotropic-3d: $(TARGET11)
	./$(TARGET11)

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
	./$(TARGET10)
	./$(TARGET11)

# Run memory bandwidth test
run-memory-bw: $(MEMORY_BW_TEST)
	./$(MEMORY_BW_TEST)

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
convolution-anisotropic-2d: $(TARGET10)
convolution-anisotropic-3d: $(TARGET11)
memory-bw: $(MEMORY_BW_TEST)

.PHONY: all clean debug-fftw run-forward run-backward run-convolution run-forward-2d run-backward-2d run-convolution-2d run-forward-3d run-backward-3d run-convolution-3d run-all run-memory-bw forward backward convolution forward-2d backward-2d convolution-2d forward-3d backward-3d convolution-3d convolution-anisotropic-2d convolution-anisotropic-3d memory-bw