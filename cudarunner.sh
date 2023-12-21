#!/bin/bash
#SBATCH --partition=cpar
#SBATCH --constraint=k20
make clean
make
time make runcuda
