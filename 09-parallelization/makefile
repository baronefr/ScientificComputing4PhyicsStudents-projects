# Compiler and flags
CC = mpicc
CFLAGS = -Xpreprocessor -fopenmp -I/opt/homebrew/Cellar/libomp/21.1.1/include -L/opt/homebrew/Cellar/libomp/21.1.1/lib -lomp -lm

# Target executable
TARGET = main

# Source files
SRC = main.c

# Default target
all: $(TARGET)

# Build the executable
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

# Build with FILE_OUTPUT defined
file: $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

# Clean up generated files
clean:
	rm -f $(TARGET)

# Run the program with example arguments
run: $(TARGET)
	mpirun -np 4 ./$(TARGET) 1000000

.PHONY: all clean run
