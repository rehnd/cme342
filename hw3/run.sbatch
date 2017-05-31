#!/bin/bash 

#SBATCH --job-name=long-run-16
#SBATCH --output=openmp.out
#SBATCH --error=openmp.err
#SBATCH --time=0:10:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=16000
#SBATCH --ntasks-per-node=16

echo "serial"
time ./serial -i input.nml

echo "hw3_1"
time ./hw3_1 -i input.nml

#echo "hw3_2"
#time ./hw3_2 -i input.nml
