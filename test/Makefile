# TempestBase directory
TEMPESTBASEDIR= ..

# Library name
LIBNAME= libhardcoreatm

# Compile with BLAS libraries
USEBLAS= true

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

# Library files to include
CFLAGS+= -I$(TEMPESTBASEDIR)/include

LDFILES= -ltempestbase -lhardcoreatm -lnetcdf -lnetcdf_c++ -lX11 -lpetsc

LDFLAGS+= -L$(TEMPESTBASEDIR)/lib -L$(NETCDF_LIBDIR) -L$(PETSC_LIBDIR) -L$(X11_LIBDIR)

##
## Build instructions
##
all: atm SWTest2 BarotropicInstabilityTest BaroclinicWaveJWTest BaroclinicWaveUMJSTest

atm:
	cd $(TEMPESTBASEDIR)/src/base; make
	cd $(TEMPESTBASEDIR)/src/atm; make

##
## Individual test case build instructions
##
SWTest2: $(BUILDDIR)/SWTest2.o $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/SWTest2.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BarotropicInstabilityTest: $(BUILDDIR)/BarotropicInstabilityTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BarotropicInstabilityTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BaroclinicWaveJWTest: $(BUILDDIR)/BaroclinicWaveJWTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaroclinicWaveJWTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BaroclinicWaveUMJSTest: $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f SWTest2
	rm -f BarotropicInstabilityTest
	rm -f BaroclinicWaveJWTest
	rm -f BaroclinicWaveUMJSTest
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

