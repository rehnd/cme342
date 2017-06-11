from pylab import *


names = ['MPI_Send\nMPI_Recv',
         'MPI_SendRecv',
         'MPI_ISend\nMPI_IRecv',
         'Coarrays',
         'Coarrays with\nDerived Data types']
x = arange(len(names))

nlines = array([25, 17, 41, 10, 9])

fs = 18

f = figure(figsize=(12,10))
bar(x,nlines)
xticks(x,names)
ylabel("Number of lines of code", fontsize=fs)
xlabel("MPI/Coarray Problem", fontsize=fs)
title("Comparison of code complexity",fontsize=fs)
show()
