# Compiler and flags
CC = gcc
CFLAGS = -O3 -ffast-math -funroll-loops -Wall -Wextra
LDFLAGS = -lfftw3 -lm

# Common source files (without the main files)
COMMON_SRC = twiddle_factors.c butterflies.c index_reversals.c fft.c convolution.c misc.c
COMMON_OBJ = $(COMMON_SRC:.c=.o)
DEPS = fft.h

# Executables
TARGET1 = bench1d_forward_fft
TARGET2 = bench1d_backward_fft
TARGET3 = bench1d_convolution
TARGETS = $(TARGET1) $(TARGET2) $(TARGET3)

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

# Compile .c to .o
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(COMMON_OBJ) bench1d_forward_fft.o bench1d_backward_fft.o bench1d_convolution.o $(TARGETS)

# Run benchmarks
run-forward: $(TARGET1)
	./$(TARGET1)

run-backward: $(TARGET2)
	./$(TARGET2)

run-convolution: $(TARGET3)
	./$(TARGET3)

run-all: $(TARGETS)
	./$(TARGET1)
	./$(TARGET2)
	./$(TARGET3)

# Individual targets (useful for building just one)
forward: $(TARGET1)
backward: $(TARGET2)
convolution: $(TARGET3)

.PHONY: all clean run-forward run-backward run-convolution run-all forward backward convolution