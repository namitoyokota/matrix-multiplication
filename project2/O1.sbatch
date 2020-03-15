#!/bin/bash -l

#SBATCH -q defg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -J O1
#SBATCH -o O1.out
#SBATCJ -e O1.%j.error

./O1