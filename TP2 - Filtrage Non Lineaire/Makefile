CC=gcc
CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: median


median: median.o pgm.o fft.o divers.o
	$(CC)  $^  -o $@ $(LDFLAGS)

%.o:%.c
	$(CC) -lm  $(CFLAGS) $<



