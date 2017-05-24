#!/bin/bash

#SBATCH -p evanreed
#SBATCH --job-name=parmetis
#SBATCH --output=metis.out
#SBATCH --error=metis.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=20

./serial 2 ../meshes/oneram6.conn ../meshes/oneram6.xyz
./serial 4 ../meshes/oneram6.conn ../meshes/oneram6.xyz
./serial 8 ../meshes/oneram6.conn ../meshes/oneram6.xyz
./serial 16 ../meshes/oneram6.conn ../meshes/oneram6.xyz
./serial 32 ../meshes/oneram6.conn ../meshes/oneram6.xyz
./serial 64 ../meshes/oneram6.conn ../meshes/oneram6.xyz
./serial 128 ../meshes/oneram6.conn ../meshes/oneram6.xyz
echo ''
./serial 2 ../meshes/falcon.conn ../meshes/falcon.xyz
./serial 4 ../meshes/falcon.conn ../meshes/falcon.xyz
./serial 8 ../meshes/falcon.conn ../meshes/falcon.xyz
./serial 16 ../meshes/falcon.conn ../meshes/falcon.xyz
./serial 32 ../meshes/falcon.conn ../meshes/falcon.xyz
./serial 64 ../meshes/falcon.conn ../meshes/falcon.xyz
./serial 128 ../meshes/falcon.conn ../meshes/falcon.xyz
echo ''
./serial 2 ../meshes/glidebak.conn ../meshes/glidebak.xyz
./serial 4 ../meshes/glidebak.conn ../meshes/glidebak.xyz
./serial 8 ../meshes/glidebak.conn ../meshes/glidebak.xyz
./serial 16 ../meshes/glidebak.conn ../meshes/glidebak.xyz
./serial 32 ../meshes/glidebak.conn ../meshes/glidebak.xyz
./serial 64 ../meshes/glidebak.conn ../meshes/glidebak.xyz
./serial 128 ../meshes/glidebak.conn ../meshes/glidebak.xyz
echo ''

