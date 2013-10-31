# TempestBase directory
TEMPESTBASEDIR= ../..

# Library name
LIBNAME= libhardcoreatm

# Compile with BLAS libraries
USEBLAS= true

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

# Local files
FILES= CubedSphereTrans.cpp \
       GridSpacing.cpp \
       PatchBox.cpp \
       ConsolidationStatus.cpp \
       Connectivity.cpp \
       Model.cpp \
       TimestepSchemeARK4.cpp \
       HorizontalDynamicsFEM.cpp \
       VerticalDynamicsFEM.cpp \
       JacobianFreeNewtonKrylov.cpp \
       FluxReconstructionFunction.cpp \
       OutputManager.cpp \
       OutputManagerComposite.cpp \
       OutputManagerReference.cpp \
       OutputManagerChecksum.cpp \
       PhysicalConstants.cpp \
       GaussLobattoQuadrature.cpp \
       GaussQuadrature.cpp \
       Grid.cpp \
       GridPatch.cpp \
       GridLonLat.cpp \
       GridGLL.cpp \
       GridCSGLL.cpp \
       GridPatchGLL.cpp \
       GridPatchCSGLL.cpp \
       Direction.cpp \
       GridData3D.cpp \
       GridData4D.cpp \
       TestCase.cpp \
	   TimeObj.cpp

##
## Build flags
##
CFLAGS+= -I$(TEMPESTBASEDIR)/include

##
## Build instructions
##
all: $(FILES:%.cpp=$(BUILDDIR)/%.o)
	rsync --checksum *.h $(TEMPESTBASEDIR)/include
	rm -f $(TEMPESTBASEDIR)/lib/$(LIBNAME).a
	ar -cqs $(TEMPESTBASEDIR)/lib/$(LIBNAME).a build/*.o

##
## Clean
##
clean:
	rm -f $(TEMPESTBASEDIR)/lib/$(LIBNAME).a
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

