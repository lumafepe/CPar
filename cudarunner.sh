#!/bin/bash
#SBATCH --partition=cpar
#SBATCH --constraint=k20
make clean
make
module load gcc/7.2.0
module load cuda/11.3.1
nvprof --profile-child-processes make runcuda
