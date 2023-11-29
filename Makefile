CC = gcc
SRC = src/
CFLAGS =-O3 -march=native -ftree-vectorize -mavx -Wall -fno-omit-frame-pointer -lm

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm ./MD*.exe

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	./MDpar.exe < inputdata.txt
