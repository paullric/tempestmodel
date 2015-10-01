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
       ConsolidationStatus.cpp \
       Connectivity.cpp \
       Model.cpp \
	   EquationSet.cpp \
       TimestepSchemeStrang.cpp \
       TimestepSchemeARK2.cpp \
       TimestepSchemeARK3.cpp \
       TimestepSchemeARK4.cpp \
       HorizontalDynamicsFEM.cpp \
       HorizontalDynamicsDG.cpp \
       VerticalDynamicsFEM.cpp \
       JacobianFreeNewtonKrylov.cpp \
       FluxCorrectionFunction.cpp \
       LinearColumnOperator.cpp \
       LinearColumnOperatorFEM.cpp \
       OutputManager.cpp \
       OutputManagerComposite.cpp \
       OutputManagerReference.cpp \
       OutputManagerChecksum.cpp \
       PhysicalConstants.cpp \
       Grid.cpp \
       GridPatch.cpp \
       GridGLL.cpp \
       GridCSGLL.cpp \
       GridCartesianGLL.cpp \
       GridPatchGLL.cpp \
       GridPatchCSGLL.cpp \
       GridPatchCartesianGLL.cpp \
       Direction.cpp \
       TestCase.cpp \
	   WorkflowProcess.cpp \
	   HeldSuarezPhysics.cpp

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

