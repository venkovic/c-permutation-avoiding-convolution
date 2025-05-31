# Compiler and flags
CC = gcc
CFLAGS = -O3 -ffast-math -funroll-loops -Wall -Wextra
LDFLAGS = -lfftw3 -lm

# Project files
SRC = twiddle_factors.c butterflies.c index_reversals.c fft.c convolution.c misc.c bench1d.c
OBJ = $(SRC:.c=.o)
DEPS = fft.h

# Output executable
TARGET = bench1d

# Default rule
all: $(TARGET)

# Link the object files into the final executable
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile .c to .o
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJ) $(TARGET)

# Run the benchmark
run: $(TARGET)
	./$(TARGET)
