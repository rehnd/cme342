#!/bin/bash 

#SBATCH --job-name=long-run-16
#SBATCH --output=hw1.out
#SBATCH --error=hw1.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=16

mpirun -np 16 ./hw1_i