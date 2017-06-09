#!/bin/bash 

#SBATCH -p evanreed
#SBATCH --job-name=long-run-16
#SBATCH --output=project.out
#SBATCH --error=project.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=2
#SBATCH --mem=48000
#SBATCH --ntasks-per-node=20

python run.py
