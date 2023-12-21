CC = gcc
NVCC = nvcc
SRC = src/
NFLAGS = -O3 -lm -arch=sm_61
CFLAGS =-O3 -march=native -ftree-vectorize -mavx -Wall -fno-omit-frame-pointer -lm

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe MDcuda.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -fopenmp -o MDpar.exe
MDcuda.exe: $(SRC)/MDcuda.cu
	module load gcc/7.2.0;\
	module load cuda/11.3.1;\
	$(NVCC) $(NFLAGS) $(SRC)MDcuda.cu -o MDcuda.exe

clean:
	rm ./MD*.exe

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	./MDpar.exe < inputdata.txt

runcuda: MDcuda.exe
	./MDcuda.exe < inputdata.txt


PROF_FLAGS = -pg

prof: $(SRC)/MDcuda.cu
	$(NVCC) $(NFLAGS) $(PROF_FLAGS) $(SRC)MDcuda.cu -lm -o prof_md

run-prof: prof
	./prof_md < inputdata.txt

graph-prof: run-prof
	gprof prof_md > main.gprof
	gprof2dot -o output.dot main.gprof
	rm gmon.out
	dot -Tpng -o output.png output.dot