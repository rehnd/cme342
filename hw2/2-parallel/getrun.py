from pylab import *


filenames = ["oneram6", "falcon", "glidebak"]
np = array([2**n for n in range(1,7)])
meshdir = "../meshes/"
runstring = """#!/bin/bash

#SBATCH -p evanreed
#SBATCH --job-name=parmetis
#SBATCH --output=metis.out
#SBATCH --error=metis.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=4
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=20

"""

for f in filenames:
    f1 = meshdir + f + ".conn"
    f2 = meshdir + f + ".xyz"
    for n in np:
        runstring += "mpirun -np " + str(n) + " ./parallel " + str(n) + " " + f1 + " " + f2 + "\n"
    runstring += "echo ''\n"

print(runstring)
