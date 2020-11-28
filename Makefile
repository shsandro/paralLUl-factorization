CC=gcc
UPCC=upcc

all:
	gcc sequential/lu_seq.c -o lu_seq.out
	upcc -network=smp parallel/lu_parallel.upc -o lu_par.out

clean:
	rm -rf *.out