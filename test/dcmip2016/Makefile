# TempestBase directory
TEMPESTBASEDIR= ../..


# Compile with BLAS libraries
USEBLAS= true

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

##
## Build instructions
##
all: atm TropicalCycloneTest BaroclinicWaveUMJSTest
 
atm:
	cd $(TEMPESTBASEDIR)/src/base; make
	cd $(TEMPESTBASEDIR)/src/atm; make
	cd interface; make 

INTERFACE_FILES= interface/main.o interface/tropical_cyclone_test.o interface/Terminator.o interface/test_baroclinic.o

##
## Individual test case build instructions
##
BaroclinicWaveUMJSTest : $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(INTERFACE_FILES) $(LDFILES) -lgfortran

TropicalCycloneTest : $(BUILDDIR)/TropicalCycloneTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/TropicalCycloneTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(INTERFACE_FILES) $(LDFILES) -lgfortran

##
## Clean
##
clean:
	rm -f BaroclinicWaveUMJSTest 
	rm -f TropicalCycloneTest 	
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

