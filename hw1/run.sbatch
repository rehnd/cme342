#!/bin/bash 

#SBATCH --job-name=long-run-16
#SBATCH --output=MPI_cm342.out
#SBATCH --error=MPI_cm342.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=16000
#SBATCH --ntasks-per-node=16

python run.py