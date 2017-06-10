import os
from pylab import *

def gettimes(probs,nps):
    times = zeros([len(probs),len(nps)])


    j = 0
    for np in nps:
        k = 0
        for prob in probs:
            fname = 'tmpdir/'+prob+'_np'+str(np)
            tmptime = genfromtxt(fname)
            times[k,j] = average(tmptime)
            k += 1
        j += 1

    return times
    


def doplots():
    fs = 22
    lw = 5
    leg = ["MPI_Send, MPI_Recv",
           "MPI_SendRecv",
           "MPI_ISend, MPI_IRecv",
           "Coarrays",
           "Coarrays w/ derived datatypes"]
    
    f = figure(figsize=(12,10))
    for i in range(len(probs)):
        plot(nps[:],times[i,:],linewidth=lw)

    xticks([1,2,4])
    legend(leg,fontsize=fs)
    xlabel('Number of processors',fontsize=fs)
    ylabel('Wall time (sec)',fontsize=fs)
    show()

    f = figure(figsize=(12,10))
    for i in range(len(probs)):
        plot(nps[:],times[i,0]/times[i,:],linewidth=lw)

    plot(nps[:],nps[:],'k--',linewidth=lw)
    xticks([1,2,4])
    legend(probs + ['ideal'], loc=2,fontsize=fs)
    xlabel('Number of processors',fontsize=fs)
    ylabel('Speedup',fontsize=fs)
    show()

    f = figure(figsize=(12,10))
    for i in range(len(probs)):
        plot(nps[:],times[i,0]/times[i,:]/nps,linewidth=lw)

    plot(nps[:],ones(len(nps)),'k--',linewidth=lw)
    ylim(0,1.3)
    xticks([1,2,4])
    legend(probs+['ideal'], loc=3,fontsize=fs)
    xlabel('Number of processors',fontsize=fs)
    ylabel('Efficiency',fontsize=fs)
    show()

if __name__ == '__main__':
    probs = ['hw1_i','hw1_ii','hw1_iii', 'coarray', 'border']
    nps = array([1,2,4])

    if ('times.dat' in os.listdir('.')):
        print("Reading times from existing times.dat")
        times = genfromtxt('times.dat')
        doplots()
    else:
        print("Generating times.dat from scratch")
        times = gettimes(probs, nps)
        savetxt('times.dat', times)
        doplots()
