# This is an commentary line in a makefile
# Start of the makefile
 
all: Terminator.o test_baroclinic.o tropical_cyclone_test.o main.o input.o

Terminator.o: Terminator.f90
	gfortran -c Terminator.f90
test_baroclinic.o: test_baroclinic.f90
	gfortran -c test_baroclinic.f90
tropical_cyclone_test.o: tropical_cyclone_test.f90
	gfortran -c tropical_cyclone_test.f90

main.o: test_baroclinic.mod Terminator.mod tropical_cyclone_test.mod main.f90
	gfortran -c main.f90
input.o: input.f90
	gfortran -c input.f90
clean:
	rm Terminator.mod Terminator.o test_baroclinic.mod test_baroclinic.o tropical_cyclone_test.o tropical_cyclone_test.mod input.o main.o test.exe
# End of the makefile
