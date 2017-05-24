# The following is a Makefile which demonstrates how to link against Metis; I
# found documentation online very misleading and unhelpful, and had to ask
# Julien how to get the damn thing working.

# The following makes the assumption that you have installed Metis to the
# directory DIR; you can change this argument for your particular installation

# If you're not a fan of Makefiles (why wouldn't you be?), you can just pattern
# match the arcane symbology in the main: call below

# Cheers -- ZDR

DIR      = /home/rehnd/parmetis-4.0.3
filename = parallel
CXX      = icpc
MPICC    = mpicxx
NPROC    = 1
LINKPARMETIS  = ${DIR}/build/Linux-x86_64/libparmetis
LINKMETIS  = ${DIR}/build/Linux-x86_64/libmetis

main:
	${MPICC} ${filename}.c -o ${filename} -std=c++11 -Wall -Wextra -g -gdwarf-4 -I ${DIR}/include/ \
	 -L ${LINKPARMETIS} -lparmetis -L ${LINKMETIS} -lmetis


clean:
	rm -f ${filename}.c~
	rm -f Makefile~
	rm -f ${filename}
	rm -rf ${filename}.dSYM
