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

function writeinput()
{
    echo "&input_parameters"
    echo "    n1  = "$1
    echo "    n2  = "$2
    echo "    np1 = "$3
    echo "    np2 = "$4
    echo "    niter = "$5
    echo "/"
    echo " "
}

mkdir tmpdir

writeinput 1024 1024 1 1 100 > input.nml
for i in {1..10}
do
    cafrun -np 1 ./hw1_i input.nml | grep Wall | awk '{print $4}' >> tmpdir/hw1_i_np1
    cafrun -np 1 ./coarray input.nml | grep Wall | awk '{print $4}' >> tmpdir/coarray_np1
    cafrun -np 1 ./border input.nml | grep Wall | awk '{print $4}' >> tmpdir/border_np1
done

writeinput 1024 1024 2 1 100 > input.nml
for i in {1..10}
do
    cafrun -np 2 ./hw1_i input.nml | grep Wall | awk '{print $4}' >> tmpdir/hw1_i_np2
    cafrun -np 2 ./coarray input.nml | grep Wall | awk '{print $4}' >> tmpdir/coarray_np2
    cafrun -np 2 ./border input.nml | grep Wall | awk '{print $4}' >> tmpdir/border_np2
done

writeinput 1024 1024 2 2 100 > input.nml
for i in {1..10}
do
    cafrun -np 4 ./hw1_i input.nml | grep Wall | awk '{print $4}' >> tmpdir/hw1_i_np4
    cafrun -np 4 ./coarray input.nml | grep Wall | awk '{print $4}' >> tmpdir/coarray_np4
    cafrun -np 4 ./border input.nml | grep Wall | awk '{print $4}' >> tmpdir/border_np4
done

