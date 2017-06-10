from pylab import *
import os

def getheader():
    
    string = """#!/bin/bash 
    
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

"""
    return string



def writeinput(np):
    # npdict[np] returns [np1,np2]
    npdict = {1:[1,1], 2:[2,1], 4:[2,2], 8:[4,2], 16:[4,4]}
    string = '&input_parameters\n' + \
        '    n1  = 2048\n' + \
        '    n2  = 2048\n' + \
        '    np1 = ' + str(npdict[np][0]) + '\n' +\
        '    np2 = ' + str(npdict[np][1]) + '\n' + \
        '    niter = 100\n/\n'

    f = open('input.nml','w')
    f.write(string)
    f.close()

def writeinput2(np):
    npdict = {1:[1,1], 2:[2,1], 4:[2,2], 8:[4,2], 16:[4,4]}
    n1 = str(1024)
    n2 = str(1024)
    np1 = str(npdict[np][0])
    np2 = str(npdict[np][1])
    niter = str(100)

    return "writeinput %s %s %s %s %s > input.nml\n" %(n1,n2,np1,np2,niter)

def run_simulations(prob,nps,nsims_to_avg, printnorms):
    times = zeros([len(probs),len(nps)])
    if (not 'tmpdir' in os.listdir('.')):
        os.system('mkdir tmpdir')
    else:
        os.system('rm -rf tmpdir; mkdir tmpdir')

    f = open('timing.sbatch', 'w')
    f.write(getheader())

    j = 0
    for np in nps:

        f.write(writeinput2(np))
        cmd = 'for i in {1..' + str(nsims_to_avg) + '}\ndo\n'
        for prob in probs:
            fname = 'tmpdir/'+prob+'_np'+str(np)

            if printnorms:
                cmd += '    cafrun -np ' + str(np) + ' ./' + prob + ' input.nml' + \
                      " | grep norm | awk '{print $2}' >> " + fname + "\n"
            else:
                cmd += '    cafrun -np ' + str(np) + ' ./' + prob + ' input.nml' + \
                      " | grep Wall | awk '{print $4}' >> " + fname + "\n"
                
        cmd += "done\n\n"
        f.write(cmd)

                # tmptime = genfromtxt(fname)
                # times[k,j] = average(tmptime)
        j += 1

    f.close()
    return times


if __name__ == '__main__':
    # Set default parameters to get:
    printnorms = False # this option prints the norms to files, instead of times
    nsims_to_avg = 10
    #probs = ['hw1_i','hw1_ii','hw1_iii', 'coarray', 'border']
    probs = ['hw1_i', 'coarray', 'border']
    nps = array([1,2,4])

    if 'times.dat' in os.listdir('.'):
        print( "Found times.dat ; plotting from this file")
        print( "To rerun, delete times.dat before re-running")
        times = genfromtxt('times.dat')
    else:
        times = run_simulations(probs,nps,nsims_to_avg,printnorms)
        savetxt('times.dat',times)

