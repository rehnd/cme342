OS := $(shell uname)

filename = serial
FLAGS = -Wall -g -gdwarf-4

ifeq ($(OS), Darwin)
  CC       = gcc
  MPICC    = mpicc
  IMETIS   = /Users/rehnd/Desktop/metis/metis-5.1.0/include
  LMETIS   = /Users/rehnd/Desktop/metis/metis-5.1.0/build/Darwin-x86_64/libmetis
else
  CC       = g++ # icpc
  MPICC    = mpicxx
  IMETIS   = /home/users/rehnd/metis-5.1.0/include
  LMETIS   = /home/users/rehnd/metis-5.1.0/build/Linux-x86_64/libmetis
endif

main:
	${CC} ${filename}.c -o ${filename} ${FLAGS} \
	  -I ${IMETIS} \
	  -L ${LMETIS} -lmetis

clean:
	rm -f ${filename}.c~
	rm -f Makefile~
	rm -f ${filename}
	rm -rf ${filename}.dSYM
