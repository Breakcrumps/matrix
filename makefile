all: src/algo.c src/io.c matrix.c
	gcc -O2 src/algo.c src/io.c matrix.c -o matrix
