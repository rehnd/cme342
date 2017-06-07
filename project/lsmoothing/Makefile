.PHONY: clean

default:
	@echo 'Usage:'
	@echo '      make <option>'
	@echo '         <option> = coarray or serial'

all: serial coarray border

serial: serial.f90
	caf -o serial serial.f90

coarray: coarray.f90
	caf -o coarray coarray.f90

border: coarray_borders.f90
	caf -o border coarray_borders.f90

clean:
	rm serial coarray border
