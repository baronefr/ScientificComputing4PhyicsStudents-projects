# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -I/opt/homebrew/Cellar/gsl/2.8/include -L/opt/homebrew/Cellar/gsl/2.8/lib -lm -lgsl -lgslcblas

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
	./$(TARGET)

.PHONY: all clean run
