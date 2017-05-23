# The following is a Makefile which demonstrates how to link against Metis; I
# found documentation online very misleading and unhelpful, and had to ask
# Julien how to get the damn thing working.

# The following makes the assumption that you have installed Metis to the
# directory DIR; you can change this argument for your particular installation

# If you're not a fan of Makefiles (why wouldn't you be?), you can just pattern
# match the arcane symbology in the main: call below

# Cheers -- ZDR

DIR      = /home/rehnd/metis-5.1.0
filename = serial
CXX      = g++ # icpc
MPICC    = mpicxx
NPROC    = 1

main:
	${CXX} ${filename}.c -o ${filename} -std=c++11 -Wall -g -gdwarf-4 -I ${DIR}/include/ \
	 -L ${DIR}/build/Linux-x86_64/libmetis/ -lmetis	


clean:
	rm -f ${filename}.c~
	rm -f Makefile~
	rm -f ${filename}
	rm -rf ${filename}.dSYM
