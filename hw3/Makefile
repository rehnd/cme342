.PHONY: clean

default:
	@echo 'Usage:'
	@echo '      make <option>'
	@echo '         <option> = 1, 2, or 3 corresponding to the question number'

all: serial 1

serial: hw3_serial.f90
	mpif90 -o serial -qopenmp hw3_serial.f90

1: hw3_1_pardo.f90
	mpif90 -o hw3_1 -qopenmp hw3_1_pardo.f90

clean:
	rm serial hw3_1 *~ *.o
