from pylab import *
import os

def writeinput(np):
    # npdict[np] returns [np1,np2]
    npdict = {1:[1,1], 2:[2,1], 4:[2,2], 8:[4,2], 16:[4,4]}
    
    string = '&input_parameters\n' + \
        '    n1  = 1024\n' + \
        '    n2  = 1024\n' + \
        '    np1 = ' + str(npdict[np][0]) + '\n' +\
        '    np2 = ' + str(npdict[np][1]) + '\n' + \
        '    niter = 100\n' + \
        '    nthreads = ' + str(np) + '\n/'
    
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
            fname = 'tmpdir/p'+prob+'_np'+str(np)
            for i in range(nsims_to_avg):
                if printnorms:
                    cmd = './hw3_' + prob + ' -i input.nml' + \
                        " | grep norm | awk '{print $3}' >> " + fname
                else:
                    cmd = './hw3_' + prob + ' -i input.nml' + \
                        " | grep time | awk '{print $3}' >> " + fname
                #print(cmd)
                os.system(cmd)


            tmptime = genfromtxt(fname)
            times[k,j] = average(tmptime[1:-1])
            k += 1
        j += 1

    return times

if __name__ == '__main__':
    # Set default parameters to get:
    printnorms = False # this option prints the norms to files, instead of times
    nsims_to_avg = 100
    probs = ['1','2']
    nps = array([1,2,4,8,16])

    if 'times.dat' in os.listdir('.'):
        print( "Found times.dat ; plotting from this file")
        print( "To rerun, delete times.dat before re-running")
        times = genfromtxt('times.dat')
    else:
        times = run_simulations(probs,nps,nsims_to_avg,printnorms)
        savetxt('times.dat',times)



    f = figure()
    plot(nps[:],times[0,:])
    plot(nps[:],times[1,:])
    legend(["Problem 1", "Problem 2"])
    xlabel('Number of processors')
    ylabel('Wall time (sec)')
    title("Homework 3 wall times")
    show()

    f = figure()
    plot(nps[:],times[0,0]/times[0,:])
    plot(nps[:],times[0,0]/times[1,:])
    plot(nps[:],nps[:],'k--')
    legend(['Problem 1','Problem 2','ideal'], loc=2)
    xlabel('Number of processors')
    ylabel('Speedup')
    show()

    f = figure()
    plot(nps[:],times[0,0]/times[0,:]/nps)
    plot(nps[:],times[1,0]/times[1,:]/nps)
    plot(nps[:],ones(len(nps)),'k--')
    ylim(0,1.3)
    legend(["Problem 1","Problem 2"], loc=4)
    xlabel('Number of processors')
    ylabel('Efficiency')
    show()
    
