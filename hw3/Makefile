.PHONY: clean

FC=ifort        # gfortran   or   ifort only
OMPLIB=qopenmp  # fopenmp for gfortran;  qopenmp for ifort

default:
	@echo 'Usage:'
	@echo '      make <option>'
	@echo '         <option> = 1, 2, or 3 corresponding to the question number'

all: serial 1 2

serial: hw3_0_serial.f90
	$(FC) -o serial -$(OMPLIB) hw3_0_serial.f90

1: hw3_1_pardo.f90
	$(FC) -o hw3_1 -$(OMPLIB) hw3_1_pardo.f90

2: hw3_2_parregions.f90
	$(FC) -o hw3_2 -$(OMPLIB) hw3_2_parregions.f90

clean:
	rm serial hw3_1 hw3_2 *~ *.o
