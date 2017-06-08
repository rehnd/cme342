.PHONY: clean

default:
	@echo 'Usage:'
	@echo '      make <option>'
	@echo '         <option> = coarray or serial'

all: serial coarray border i ii iii

serial: serial.f90
	caf -o serial serial.f90

coarray: coarray.f90
	caf -o coarray coarray.f90

border: coarray_borders.f90
	caf -o border coarray_borders.f90

i: hw1_1_SEND_RECV.f90
	mpif90 -o hw1_i hw1_1_SEND_RECV.f90

ii: hw1_2_SENDRECV.f90
	mpif90 -o hw1_ii hw1_2_SENDRECV.f90

iii: hw1_3_ISEND_IRECV.f90
	mpif90 -o hw1_iii hw1_3_ISEND_IRECV.f90

clean:
	rm serial coarray border hw1_i hw1_ii hw1_iii

