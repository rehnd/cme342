.PHONY: clean

FLAGS = -O3
default:
	@echo 'Usage:'
	@echo '      make <option>'
	@echo '         <option> = coarray or serial'

all: serial coarray border i ii iii

serial: serial.f90
	gfortran -o serial $(FLAGS) serial.f90

coarray: coarray.f90
	caf -o coarray $(FLAGS) coarray.f90

border: coarray_borders.f90
	caf -o border $(FLAGS) coarray_borders.f90

i: hw1_1_SEND_RECV.f90
	caf -o hw1_i $(FLAGS) hw1_1_SEND_RECV.f90

ii: hw1_2_SENDRECV.f90
	caf -o hw1_ii $(FLAGS) hw1_2_SENDRECV.f90

iii: hw1_3_ISEND_IRECV.f90
	caf -o hw1_iii $(FLAGS) hw1_3_ISEND_IRECV.f90

clean:
	rm serial coarray border hw1_i hw1_ii hw1_iii

