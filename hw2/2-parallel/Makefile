OS := $(shell uname)

filename = parallel
FLAGS = -Wall -Wextra -g -gdwarf-4 -lm

CC        = gcc
MPICC     = mpicc

ifeq ($(OS),Darwin)
  DIR       = /Users/rehnd/Desktop/metis
  IMETIS    = ${DIR}/metis-5.1.0/include
  IPMETIS   = ${DIR}/parmetis-4.0.3/include
  LMETIS    = ${DIR}/metis-5.1.0/build/Darwin-x86_64/libmetis
  LPMETIS   = ${DIR}/parmetis-4.0.3/build/Darwin-x86_64/libparmetis
else # assume OS == Linux
  DIR      = /home/users/rehnd
  IMETIS   = ${DIR}/metis-5.1.0/include
  IPMETIS  = ${DIR}/parmetis-4.0.3/include
  LMETIS   = ${DIR}/metis-5.1.0/build/Linux-x86_64/libmetis
  LPMETIS  = ${DIR}/parmetis-4.0.3/build/Linux-x86_64/libparmetis
endif

main:
	${MPICC} ${filename}.c -o ${filename} ${FLAGS} \
	 -I ${IPMETIS} \
	 -L ${LPMETIS} -lparmetis \
	 -I ${IMETIS} \
	 -L ${LMETIS} -lmetis

clean:
	rm -f ${filename}.c~
	rm -f Makefile~
	rm -f ${filename}
	rm -rf ${filename}.dSYM
