# Tempest base directory
TEMPESTBASEDIR= ../../base

# TempestBase directory
BASE_INCLUDEDIR= $(TEMPESTBASEDIR)/include
BASE_LIBDIR= $(TEMPESTBASEDIR)/lib

# PetSc include directory
PETSC_INCLUDEDIR= /usr/local/include
PETSC_LIBDIR= /usr/local/lib

# NETCDF directories
NETCDF_INCLUDEDIR=/usr/local/include
NETCDF_LIBDIR=/usr/local/lib

# X11 lib directory
X11_LIBDIR=/usr/X11/lib

# Use BLAS
USEBLAS= true

# Library files to include
LDFILES= -ltempestbase -lnetcdf -lnetcdf_c++ -lX11 -lpetsc

# Local files
FILES= CubedSphereTrans.cpp GridSpacing.cpp PatchBox.cpp \
       ConsolidationStatus.cpp Connectivity.cpp \
       Model.cpp \
       TimestepSchemeARK4.cpp \
       HorizontalDynamicsFEM.cpp \
       VerticalDynamicsFEM.cpp \
       JacobianFreeNewtonKrylov.cpp \
       OutputManager.cpp OutputManagerComposite.cpp \
       OutputManagerReference.cpp OutputManagerChecksum.cpp \
       PhysicalConstants.cpp \
	   GaussLobattoQuadrature.cpp GaussQuadrature.cpp \
       Grid.cpp GridPatch.cpp GridLonLat.cpp GridCSGLL.cpp \
       Direction.cpp GridData3D.cpp GridData4D.cpp TestCase.cpp

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

CFLAGS+= -I$(BASE_INCLUDEDIR) -I$(NETCDF_INCLUDEDIR) -I$(PETSC_INCLUDEDIR) -DNETCDFENABLED
LDFLAGS+= -L$(BASE_LIBDIR) -L$(NETCDF_LIBDIR) -L$(PETSC_LIBDIR) -L$(X11_LIBDIR)

##
## Build instructions
##
all: SWTest2 BarotropicInstabilityTest BaroclinicWaveJWTest BaroclinicWaveUMJSTest

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
	rm -f SWTest2 *.o
	rm -f BaroclinicWaveJWTest *.o
	rm -f BaroclinicWaveUMJSTest *.o
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

