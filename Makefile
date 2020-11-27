CC=gcc
UPCC=upcc

all:
	gcc sequential/lu_seq.c -o lu_seq.out

clean:
	rm -rf *.out