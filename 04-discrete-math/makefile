# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -lm

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
	$(CC) $(CFLAGS) -DFILE_OUTPUT=\"samples.txt\" -o $(TARGET) $(SRC)

# Clean up generated files
clean:
	rm -f $(TARGET)

# Run the program with example arguments
run: $(TARGET)
	./$(TARGET) 100000 0 1.5707963267948966

.PHONY: all clean run
