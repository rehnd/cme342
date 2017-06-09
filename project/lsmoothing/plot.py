import os
from pylab import *


def doplots():
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

if __name__ == '__main__':
    probs = ['hw1_i','hw1_ii','hw1_iii', 'coarray', 'border']
    nps = array([1,2,4,8,16])


    if ('times.dat' in os.listdir('.')):
        times = genfromtxt('times.dat')
        doplots()
