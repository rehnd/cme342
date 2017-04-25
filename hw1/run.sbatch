#!/bin/bash 

#SBATCH --job-name=long-run-16
#SBATCH --output=mpiout.out
#SBATCH --error=mpierr.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=16

mpirun -np 16 ./hw1_i