CC=gcc

INSTALL_DIR=~/bin

all: simplematch superalign

clean: 
	rm -rf *.o
	rm -rf *.dSYM
	
distclean:
	rm -rf superalign
	rm -rf *.o
	rm -rf *.dSYM

simplematch: arrays.o matchup.o sort.o error.o simplealign.o mt19937ar.o simplematch.c
	$(CC) -g -o simplematch simplematch.c arrays.o matchup.o sort.o error.o simplealign.o mt19937ar.o -lm

superalign: arrays.o matchup.o sort.o error.o simplealign.o mt19937ar.o superalign.c
	$(CC) -g -o superalign superalign.c arrays.o matchup.o sort.o error.o  simplealign.o mt19937ar.o -lm

install:
	cp superalign simplematch ${INSTALL_DIR}