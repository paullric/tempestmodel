# TempestBase directory
TEMPESTBASEDIR= ..

# Compile with BLAS libraries
USEBLAS= true

# Library files to include
CFLAGS+= -I$(TEMPESTBASEDIR)/include

TEMPESTLIBS= \
    $(TEMPESTBASEDIR)/lib/libtempestbase.a \
	$(TEMPESTBASEDIR)/lib/libhardcoreatm.a

LDFILES= -lhardcoreatm -ltempestbase -lnetcdf -lnetcdf_c++ -lX11 -lpetsc

LDFLAGS+= -L$(TEMPESTBASEDIR)/lib -L$(NETCDF_LIBDIR) -L$(PETSC_LIBDIR) -L$(X11_LIBDIR)

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

##
## Build instructions
##
all: atm SWTest2 MountainRossbyTest BarotropicInstabilityTest BaroclinicWaveJWTest BaroclinicWaveUMJSTest MountainWaveSphereTest InertiaGravityWaveTest InertialGravityCartesianXZTest ThermalBubbleCartesianTest ShallowWaterEddyTest ScharMountainCartesianTest HydrostaticMountainCartesianTest NonHydroMountainCartesianTest

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

BarotropicInstabilityTest: $(BUILDDIR)/BarotropicInstabilityTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BarotropicInstabilityTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

ShallowWaterEddyTest: $(BUILDDIR)/ShallowWaterEddyTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/ShallowWaterEddyTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BaroclinicWaveJWTest: $(BUILDDIR)/BaroclinicWaveJWTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaroclinicWaveJWTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BaroclinicWaveUMJSTest: $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

MountainWaveSphereTest: $(BUILDDIR)/MountainWaveSphereTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/MountainWaveSphereTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

InertiaGravityWaveTest: $(BUILDDIR)/InertiaGravityWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/InertiaGravityWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

ThermalBubbleCartesianTest: $(BUILDDIR)/ThermalBubbleCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/ThermalBubbleCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

InertialGravityCartesianXZTest: $(BUILDDIR)/InertialGravityCartesianXZTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/InertialGravityCartesianXZTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

ScharMountainCartesianTest: $(BUILDDIR)/ScharMountainCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/ScharMountainCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

HydrostaticMountainCartesianTest: $(BUILDDIR)/HydrostaticMountainCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/HydrostaticMountainCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

NonHydroMountainCartesianTest: $(BUILDDIR)/NonHydroMountainCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/NonHydroMountainCartesianTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f SWTest2
	rm -f BarotropicInstabilityTest
	rm -f ShallowWaterEddyTest
	rm -f BaroclinicWaveJWTest
	rm -f BaroclinicWaveUMJSTest
	rm -f IntertialGravityCartesianXZTest
	rm -f ScharMountainCartesianTest
	rm -f HydrostaticMountainCartesianTest
	rm -f NonHydroMountainCartesianTest
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

