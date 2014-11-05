# TempestBase directory
TEMPESTBASEDIR= ../..

# Compile with BLAS libraries
USEBLAS= true

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

##
## Build instructions
##
all: atm SWTest2 MountainRossbyTest RossbyHaurwitzWaveTest BarotropicInstabilityTest

atm:
	cd $(TEMPESTBASEDIR)/src/base; make
	cd $(TEMPESTBASEDIR)/src/atm; make

##
## Individual test case build instructions
##
SWTest2: $(BUILDDIR)/SWTest2.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/SWTest2.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

MountainRossbyTest: $(BUILDDIR)/MountainRossbyTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/MountainRossbyTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

RossbyHaurwitzWaveTest: $(BUILDDIR)/RossbyHaurwitzWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/RossbyHaurwitzWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BarotropicInstabilityTest: $(BUILDDIR)/BarotropicInstabilityTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BarotropicInstabilityTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f SWTest2
	rm -f MountainRossbyTest
	rm -f RossbyHaurwitzWaveTest
	rm -f BarotropicInstabilityTest
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

