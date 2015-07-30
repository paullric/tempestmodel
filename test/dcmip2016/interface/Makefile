# This is an commentary line in a makefile
# Start of the makefile
 
all:  baroclinic_wave_test.o tropical_cyclone_test.o mesoscale_storm_test.o kessler.o

#Terminator.o: Terminator.f90
#	gfortran -c Terminator.f90
baroclinic_wave_test.o: baroclinic_wave_test.f90
	gfortran -c baroclinic_wave_test.f90
tropical_cyclone_test.o: tropical_cyclone_test.f90
	gfortran -c tropical_cyclone_test.f90
mesoscale_storm_test.o: mesoscale_storm_test.f90
	gfortran -c mesoscale_storm_test.f90
kessler.o: kessler.f90
	gfortran -c kessler.f90

#main.o: test_baroclinic.mod tropical_cyclone_test.mod mesoscale_storm_test.mod kessler.o main.f90
#	gfortran -c main.f90
	

clean:
	rm *.o *.mod
# End of the makefile
