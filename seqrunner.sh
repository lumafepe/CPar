#!/bin/bash
#SBATCH --partition=cpar
#SBATCH --cpus-per-task=1
make clean
make
perf stat -r 3 -M cpi,instructions -e branch-misses,L1-dcache-load-misses,cycles,duration_time,mem-loads,mem-stores make runseq
# time make runseq