#!/bin/bash -l

#SBATCH -q defg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -J O3
#SBATCH -o O3.out
#SBATCJ -e O3.%j.error

./O3