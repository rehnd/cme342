#!/bin/bash

#SBATCH -p evanreed
#SBATCH --job-name=parmetis
#SBATCH --output=metis.out
#SBATCH --error=metis.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=4
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=20

mpirun -np 2 ./parallel 2 ../meshes/oneram6.conn ../meshes/oneram6.xyz
mpirun -np 4 ./parallel 4 ../meshes/oneram6.conn ../meshes/oneram6.xyz
mpirun -np 8 ./parallel 8 ../meshes/oneram6.conn ../meshes/oneram6.xyz
mpirun -np 16 ./parallel 16 ../meshes/oneram6.conn ../meshes/oneram6.xyz
mpirun -np 32 ./parallel 32 ../meshes/oneram6.conn ../meshes/oneram6.xyz
echo ''
mpirun -np 2 ./parallel 2 ../meshes/falcon.conn ../meshes/falcon.xyz
mpirun -np 4 ./parallel 4 ../meshes/falcon.conn ../meshes/falcon.xyz
mpirun -np 8 ./parallel 8 ../meshes/falcon.conn ../meshes/falcon.xyz
mpirun -np 16 ./parallel 16 ../meshes/falcon.conn ../meshes/falcon.xyz
mpirun -np 32 ./parallel 32 ../meshes/falcon.conn ../meshes/falcon.xyz
echo ''
mpirun -np 2 ./parallel 2 ../meshes/glidebak.conn ../meshes/glidebak.xyz
mpirun -np 4 ./parallel 4 ../meshes/glidebak.conn ../meshes/glidebak.xyz
mpirun -np 8 ./parallel 8 ../meshes/glidebak.conn ../meshes/glidebak.xyz
mpirun -np 16 ./parallel 16 ../meshes/glidebak.conn ../meshes/glidebak.xyz
mpirun -np 32 ./parallel 32 ../meshes/glidebak.conn ../meshes/glidebak.xyz
echo ''

