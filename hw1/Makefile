.PHONY: clean

default:
	@echo 'Usage:'
	@echo '      make <option>'
	@echo '         <option> = i, ii, iii, or serial, corresponding to the question number'

all: serial i ii iii

serial: hw1_0_serial.f90
	mpif90 -o serial hw1_0_serial.f90

i: hw1_1_SEND_RECV.f90
	mpif90 -o hw1_i hw1_1_SEND_RECV.f90

ii: hw1_2_SENDRECV.f90
	mpif90 -o hw1_ii hw1_2_SENDRECV.f90

iii: hw1_3_ISEND_IRECV.f90
	mpif90 -o hw1_iii hw1_3_ISEND_IRECV.f90

clean:
	rm serial hw1_i hw1_ii hw1_iii