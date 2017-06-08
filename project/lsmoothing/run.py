from pylab import *
import os

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

def run_simulations(prob,nps,nsims_to_avg, printnorms):
    times = zeros([len(probs),len(nps)])
    if (not 'tmpdir' in os.listdir('.')):
        os.system('mkdir tmpdir')
    else:
        os.system('rm -rf tmpdir; mkdir tmpdir')

    j = 0
    for np in nps:
        k = 0
        writeinput(np)
        for prob in probs:
            fname = 'tmpdir/'+prob+'_np'+str(np)
            for i in range(nsims_to_avg):
                if printnorms:
                    cmd = 'cafrun -np ' + str(np) + ' ./' + prob + ' input.nml' + \
                        " | grep norm | awk '{print $2}' >> " + fname
                else:
                    cmd = 'cafrun -np ' + str(np) + ' ./' + prob + ' input.nml' + \
                        " | grep Wall | awk '{print $4}' >> " + fname
                os.system(cmd)

            tmptime = genfromtxt(fname)
            times[k,j] = average(tmptime)
            k += 1
        j += 1

    return times

if __name__ == '__main__':
    # Set default parameters to get:
    printnorms = False # this option prints the norms to files, instead of times
    nsims_to_avg = 10
    probs = ['hw1_i','hw1_ii','hw1_iii', 'coarray', 'border']
    nps = array([1,2,4])

    if 'times.dat' in os.listdir('.'):
        print "Found times.dat ; plotting from this file"
        print "To rerun, delete times.dat before re-running"
        times = genfromtxt('times.dat')
    else:
        times = run_simulations(probs,nps,nsims_to_avg,printnorms)
        savetxt('times.dat',times)



    f = figure()
    for i in range(len(probs)):
        plot(nps[:],times[i,:])

    legend(probs)
    xlabel('Number of processors')
    ylabel('Wall time (sec)')
    show()

    f = figure()
    for i in range(len(probs)):
        plot(nps[:],times[i,0]/times[i,:])

    plot(nps[:],nps[:],'k--')
    legend(probs + ['ideal'], loc=2)
    xlabel('Number of processors')
    ylabel('Speedup')
    show()

    f = figure()
    for i in range(len(probs)):
        plot(nps[:],times[i,0]/times[i,:]/nps)

    plot(nps[:],ones(len(nps)),'k--')
    ylim(0,1.3)
    legend(probs+['ideal'], loc=4)
    xlabel('Number of processors')
    ylabel('Efficiency')
    show()
    
