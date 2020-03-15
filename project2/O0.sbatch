#!/bin/bash -l

#SBATCH -q defg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -J O0
#SBATCH -o O0.out
#SBATCJ -e O0.%j.error

./O0