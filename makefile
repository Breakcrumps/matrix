all: src/algo.c src/io.c matrix.c
	gcc -O3 -march=x86-64-v3 src/algo.c src/io.c matrix.c -o matrix
