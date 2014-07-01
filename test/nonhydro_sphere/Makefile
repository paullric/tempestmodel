# TempestBase directory
TEMPESTBASEDIR= ../..

# Compile with BLAS libraries
USEBLAS= true

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

##
## Build instructions
##
all: atm BaroclinicWaveJWTest BaroclinicWaveUMJSTest MountainWaveSphereTest ScharMountainSphereTest InertiaGravityWaveTest BaldaufGravityWaveTest MountainRossby3DTest StationaryMountainFlowTest HeldSuarezTest
 
atm:
	cd $(TEMPESTBASEDIR)/src/base; make
	cd $(TEMPESTBASEDIR)/src/atm; make

##
## Individual test case build instructions
##
BaroclinicWaveJWTest: $(BUILDDIR)/BaroclinicWaveJWTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaroclinicWaveJWTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BaroclinicWaveUMJSTest: $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaroclinicWaveUMJSTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

MountainWaveSphereTest: $(BUILDDIR)/MountainWaveSphereTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/MountainWaveSphereTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

ScharMountainSphereTest: $(BUILDDIR)/ScharMountainSphereTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/ScharMountainSphereTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

InertiaGravityWaveTest: $(BUILDDIR)/InertiaGravityWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/InertiaGravityWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

BaldaufGravityWaveTest: $(BUILDDIR)/BaldaufGravityWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/BaldaufGravityWaveTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

MountainRossby3DTest: $(BUILDDIR)/MountainRossby3DTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/MountainRossby3DTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

StationaryMountainFlowTest: $(BUILDDIR)/StationaryMountainFlowTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/StationaryMountainFlowTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

HeldSuarezTest: $(BUILDDIR)/HeldSuarezTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/HeldSuarezTest.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f BaroclinicWaveJWTest
	rm -f BaroclinicWaveUMJSTest
	rm -f MountainWaveSphereTest
	rm -f InertiaGravityWaveTest
	rm -f BaldaufGravityWaveTest
	rm -f MountainRossby3DTest
	rm -f StationaryMountainFlowTest
	rm -f HeldSuarezTest
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

