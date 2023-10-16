CC = gcc
SRC = src/
CFLAGS = -O3 -fopt-info-vec-missed -march=native -ftree-vectorize -mavx

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm ./MD.exe

run:
	./MD.exe < inputdata.txt

# Compiling for performance testing.

PROF_FLAGS = -pg

prof: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(SRC)MD.cpp -lm -o prof_md

run-prof: prof
	./prof_md < inputdata.txt

graph-prof: run-prof
	gprof prof_md > main.gprof
	gprof2dot -o output.dot main.gprof
	rm gmon.out
	dot -Tpng -o output.png output.dot

# Compiling for manual vectorization.

VECT_FLAGS = -O3 -mavx -march=native -mtune=native

magic: $(SRC)/main.cpp
	$(CC) $(VECT_FLAGS) $(SRC)main.cpp -lm -o md_vect

cast: 
	./md_vect < inputdata.txt

vanish:
	rm cp_*
	rm md_vect

