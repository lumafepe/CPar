CC = gcc
SRC = src/
CFLAGS = -O3 -march=native -ftree-vectorize -mavx -Wall # -fopt-info-vec-all

.DEFAULT_GOAL = md.exe

md.exe: $(SRC)md.cpp
	$(CC) $(CFLAGS) $(SRC)md.cpp -lm -o MD.exe

clean:
	rm ./MD.exe

run:
	./MD.exe < inputdata.txt

# Compiling for performance testing.

PROF_FLAGS = -pg

prof: $(SRC)/md.cpp
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(SRC)md.cpp -lm -o prof_md

run-prof: prof
	./prof_md < inputdata.txt

graph-prof: run-prof
	gprof prof_md > main.gprof
	gprof2dot -o output.dot main.gprof
	rm gmon.out
	dot -Tpng -o output.png output.dot

