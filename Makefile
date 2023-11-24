CC = gcc
SRC = src
CFLAGS = -O3 -march=native -ftree-vectorize -mavx -Wall -fno-omit-frame-pointer -lm# -fopt-info-vec-all

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
	$(CC) $(CFLAGS) $< -o $@
MDpar.exe: $(SRC)/MDpar.cpp
	$(CC) $(CFLAGS) -fopenmp $< -o $@

clean:
	rm ./MD*.exe

runseq:
	./MDseq.exe < inputdata.txt

runpar:
	./MDpar.exe < inputdata.txt

# Compiling for performance testing.

PROF_FLAGS = -pg
profSeq: $(SRC)/MDseq.cpp
	$(CC) $(CFLAGS) $(PROF_FLAGS) $< -lm -o prof_md

run-profSeq: profSeq
	./prof_md < inputdata.txt

graph-profSeq: run-profSeq
	gprof prof_md > main.gprof
	gprof2dot -o output.dot main.gprof
	rm gmon.out
	dot -Tpng -o output.png output.dot

profPar: $(SRC)/MDpar.cpp
	$(CC) $(CFLAGS) $(PROF_FLAGS) $< -lm -o prof_md

run-profPar: profSeq
	./prof_md < inputdata.txt

graph-profPar: run-profSeq
	gprof prof_md > main.gprof
	gprof2dot -o output.dot main.gprof
	rm gmon.out
	dot -Tpng -o output.png output.dot

