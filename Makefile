CC=gcc

INSTALL_DIR=../../bin

all: superalign

# .o : arrays.c matchup.c sort.c error.c mt19937ar.c
#	$(CC) $(CFLAGS) -c $< -o $@
# .o.c:
# 	(CC) $(CFLAGS) $@.c -o $@

clean: 
	rm -rf *.o
	rm -rf *.dSYM
	
distclean:
	rm -rf superalign
	rm -rf *.o
	rm -rf *.dSYM

superalign: arrays.o matchup.o sort.o error.o  simplealign.o mt19937ar.o superalign.c
	$(CC) -g -o superalign superalign.c arrays.o matchup.o sort.o error.o  simplealign.o mt19937ar.o -lm

install:
	cp superalign ${INSTALL_DIR}