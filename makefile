SHELL := powershell.exe
.SHELLFLAGS := -NoProfile -Command

CC = gcc
CFLAGS = -O3 -march=x86-64-v3 -flto=auto -fno-math-errno -fno-trapping-math

TARGET = matrix

all: src/io.c src/decomps.c src/det.c src/eigen.c src/modified_mats.c src/ops.c src/rot.c src/solvers.c matrix.c
	$(CC) $(CFLAGS) -fwhole-program $^ -o $(TARGET)
sep: build/io.o build/decomps.o build/det.o build/eigen.o build/modified_mats.o build/ops.o build/rot.o build/solvers.o matrix.c
	$(CC) $(CFLAGS) -fwhole-program $^ -o $(TARGET)
build/algo.o: src/algo.c
	$(CC) -c $(CFLAGS) $< -o build/algo.o
build/decomps.o: src/decomps.c
	$(CC) -c $(CFLAGS) $< -o build/decomps.o
build/det.o: src/det.c
	$(CC) -c $(CFLAGS) $< -o build/det.o
build/eigen.o: src/eigen.c
	$(CC) -c $(CFLAGS) $< -o build/eigen.o
build/modified_mats.o: src/modified_mats.c
	$(CC) -c $(CFLAGS) $< -o build/modified_mats.o
build/ops.o: src/ops.c
	$(CC) -c $(CFLAGS) $< -o build/ops.o
build/rot.o: src/rot.c
	$(CC) -c $(CFLAGS) $< -o build/rot.o
build/solvers.o: src/solvers.c
	$(CC) -c $(CFLAGS) $< -o build/solvers.o
build/io.o: src/io.c
	$(CC) -c $(CFLAGS) $< -o build/io.o
clean:
	rm -Force build/*, *.exe
