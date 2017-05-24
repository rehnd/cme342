from pylab import *
import os

rfile = 'metis.out'

filenames = ['oneram6', 'falcon', 'glidebak']
np = array([2**n for n in range(1,8)])

if ('results' in os.listdir('.')):
    os.system("rm -rf results")
os.system('mkdir results')

for f in filenames:
    string1 = "cat " + rfile + " | grep " + f + " | grep time | awk '{print $4}' >> " + 'results/' + f + ".time"
    string2 = "cat " + rfile + " | grep " + f + " | grep edge | awk '{print $4}' >> " + 'results/' + f + ".edge"
    os.system(string1)
    os.system(string2)
    #print(string1)
    #print(string2)

times = zeros([len(np), len(filenames)])
edgecuts = zeros([len(np), len(filenames)])
for j in range(len(filenames)):
    ftname = 'results/' + filenames[j] + '.time'
    fename = 'results/' + filenames[j] + '.edge'
    t = genfromtxt(ftname)
    e = genfromtxt(fename)
    times[:,j]    = t
    edgecuts[:,j] = e

print(times)
print(edgecuts)

f = figure()
plot(np, times[:,0])
plot(np, times[:,1])
plot(np, times[:,2])
title('METIS_PartMeshDual Runtime')
xlabel('Number of Partitions')
ylabel('Time (sec)')
grid()
legend(filenames)
savefig('results/times.png',dpi=100)
show()

f = figure()
plot(np, edgecuts[:,0])
plot(np, edgecuts[:,1])
plot(np, edgecuts[:,2])
title('METIS_PartMeshDual Edgecuts (comm cost)')
xlabel('Number of Partitions')
ylabel('Time (sec)')
grid()
legend(filenames)
savefig('results/edgecuts.png',dpi=100)
show()
