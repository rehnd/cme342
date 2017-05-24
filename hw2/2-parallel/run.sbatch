#!/bin/bash 

#SBATCH --job-name=parmetis
#SBATCH --output=pmetis.out
#SBATCH --error=pmetis.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=2
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=16

mpirun -np 2 ./parallel 2 ../meshes/oneram6.conn ../meshes/oneram6.xyz
