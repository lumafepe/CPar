#!/bin/bash
module load gcc/11.2.0
make clean
make
for i in {1..40}
do
    export OMP_NUM_THREADS=$i
    echo "Threads $i"
    perf stat -r 3 -M cpi,instructions -e branch-misses,L1-dcache-load-misses,cycles,duration_time,mem-loads,mem-stores make runpar
done