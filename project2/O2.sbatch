#!/bin/bash -l

#SBATCH -q defg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -J O2
#SBATCH -o O2.out
#SBATCJ -e O2.%j.error

./O2
