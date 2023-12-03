CC=gcc
LFLAGS=-lm
CFLAGS=-Wall -Ofast -march=native -std=c99 -fopenmp

all: model

model: model.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	- /bin/rm model.o model
