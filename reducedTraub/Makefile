CC = gcc
CFLAGS = -std=c99 -Wall -O2

all: reducedTraub

reducedTraub: reducedTraub.o
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f *.o reducedTraub
