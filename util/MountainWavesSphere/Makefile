# TempestBase directory
TEMPESTBASEDIR= ../..

# Compile with BLAS libraries
USEBLAS= true

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

##
## Build instructions
##
all: atm ComputeWaveModes GenerateMountainResponse GenerateWaveTopography GenerateScharTopography

atm:
	cd $(TEMPESTBASEDIR)/src/base; make
	cd $(TEMPESTBASEDIR)/src/atm; make

##
## Individual test case build instructions
##
ComputeWaveModes: $(BUILDDIR)/ComputeWaveModes.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/ComputeWaveModes.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateMountainResponse: $(BUILDDIR)/GenerateMountainResponse.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/GenerateMountainResponse.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateWaveTopography: $(BUILDDIR)/GenerateWaveTopography.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/GenerateWaveTopography.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateScharTopography: $(BUILDDIR)/GenerateScharTopography.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/GenerateScharTopography.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

ExtractSurface: $(BUILDDIR)/ExtractSurface.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(TEMPESTLIBS)
	$(CC) $(LDFLAGS) -o $@ $(BUILDDIR)/ExtractSurface.o $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f ComputeWaveModes
	rm -f GenerateMountainResponse
	rm -f GenerateWaveTopography
	rm -f GenerateScharTopography
	rm -f ExtractSurface
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

