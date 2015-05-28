# This is an commentary line in a makefile
# Start of the makefile
 
all: Terminator.o test_baroclinic.o tropical_cyclone_test.o kessler.o main.o 

Terminator.o: Terminator.f90
	gfortran -c Terminator.f90
test_baroclinic.o: test_baroclinic.f90
	gfortran -c test_baroclinic.f90
tropical_cyclone_test.o: tropical_cyclone_test.f90
	gfortran -c tropical_cyclone_test.f90
kessler.o: kessler.f90
	gfortran -c kessler.f90

main.o: test_baroclinic.mod Terminator.mod tropical_cyclone_test.mod kessler.o main.f90
	gfortran -c main.f90
	

clean:
	rm Terminator.mod Terminator.o test_baroclinic.mod test_baroclinic.o tropical_cyclone_test.o tropical_cyclone_test.mod kessler.o main.o
# End of the makefile
