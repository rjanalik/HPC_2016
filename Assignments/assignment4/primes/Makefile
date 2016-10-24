# The following makefile works both on CUB
# Use the debug.primes executable only to debug your programs. For performance
# measurement, you should use the executable compiled with the -O flag.

CC = gcc
CFLAGS = -fopenmp -O3
DEBUGFLAGS = -fopenmp -g

all: primes primes_debug
debug: primes_debug

.PHONY: clean

clean:
	rm -rf primes primes_debug *~ *.o

primes: primes.c
	$(CC) $(CFLAGS) $< -o $@

primes_debug: primes.c
	$(CC) $(DEBUGFLAGS) $< -o $@

