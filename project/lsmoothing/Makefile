.PHONY: clean

default:
	@echo 'Usage:'
	@echo '      make <option>'
	@echo '         <option> = coarray or serial'

all: serial coarray

serial: serial.f90
	caf -o serial serial.f90

coarray: coarray.f90
	caf -o coarray coarray.f90

clean:
	rm serial coarray
