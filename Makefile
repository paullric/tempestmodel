##
## Build instructions
##
all:
	cd src/base; make
	cd src/atm; make
	cd test; make

##
## Clean
##
clean:
	cd src/base; make clean
	cd src/atm; make clean
	cd test; make clean
	rm -f include/*.h

# DO NOT DELETE
