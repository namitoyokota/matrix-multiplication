#!/bin/bash -l

#SBATCH -q defg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 00:20:00
#SBATCH -J project1
#SBATCH -o project1.out
#SBATCJ -e project1.%j.error

./project1