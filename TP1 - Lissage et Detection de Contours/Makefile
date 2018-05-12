CC=gcc
CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: essai1 essai2 convolgauss

essai1: essai1.o pgm.o
	$(CC) $^ -o $@ $(LDFLAGS)

essai2: essai2.o pgm.o fft.o
	$(CC)  $^  -o $@ $(LDFLAGS)

convolgauss: convolgauss.o pgm.o fft.o divers.o
	$(CC)  $^  -o $@ $(LDFLAGS)

%.o:%.c
	$(CC) -lm  $(CFLAGS) $<
