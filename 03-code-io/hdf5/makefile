CC = gcc
CFLAGS = -Wall -Wextra -O2 -L/usr/lib64 -lhdf5_hl -lhdf5
TARGETS = gen.x dax.x

all: $(TARGETS)

gen.x: generator.c
	$(CC) $(CFLAGS) -o gen.x generator.c

dax.x: daxpy_from_config.c
	$(CC) $(CFLAGS) -o dax.x daxpy_from_config.c

run: all
	./gen.x 10 vector
	./dax.x

clean:
	rm -f $(TARGETS) *.o

cleandata:
	rm -f *.hd5

.PHONY: all clean