CC = gcc
CFLAGS = -Wall -O2 -lgsl -lgslcblas
TARGET = gsl.axpy
SRC = gsl.axpy.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

run: all
	./$(TARGET)

clean:
	rm -f $(TARGET)
