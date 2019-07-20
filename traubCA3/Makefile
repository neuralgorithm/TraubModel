CC = gcc
CFLAGS = -std=c99 -Wall -O2

all: traubCA3

traubCA3: traubCA3.o
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f *.o traubCA3
