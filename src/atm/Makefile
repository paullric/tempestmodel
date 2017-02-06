# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTBASEDIR= ../..

# Load Makefile framework. 
include $(TEMPESTBASEDIR)/mk/framework.make

FILES= CubedSphereTrans.cpp \
       GridSpacing.cpp \
       ConsolidationStatus.cpp \
       Connectivity.cpp \
       Model.cpp \
       EquationSet.cpp \
       TimestepSchemeStrang.cpp \
       TimestepSchemeARS222.cpp \
       TimestepSchemeARS232.cpp \
       TimestepSchemeARK232.cpp \
       TimestepSchemeGARK2.cpp \
       TimestepSchemeARS343.cpp \
       TimestepSchemeARS443.cpp \
       TimestepSchemeSSP3332.cpp \
       TimestepSchemeARK4.cpp \
       TimestepSchemeSplitExp.cpp \
       HorizontalDynamicsFEM.cpp \
       HorizontalDynamicsFEMV2.cpp \
       VerticalDynamicsFEM.cpp \
       VerticalDynamicsFEMV2.cpp \
       VerticalDynamicsSchur.cpp \
       JacobianFreeNewtonKrylov.cpp \
       FluxCorrectionFunction.cpp \
       LinearColumnOperator.cpp \
       LinearColumnOperatorFEM.cpp \
       OutputManager.cpp \
       OutputManagerComposite.cpp \
       OutputManagerReference.cpp \
       OutputManagerChecksum.cpp \
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

LIBNAME= libhardcoreatm

LIB_TARGET= $(TEMPESTBASEDIR)/lib/$(LIBNAME).a

.PHONY: all clean

# Build rules.
all: $(LIB_TARGET)

$(LIB_TARGET): $(FILES:%.cpp=$(BUILDDIR)/%.o)
	rm -f $(LIB_TARGET)
	ar -cqs $(LIB_TARGET) build/*.o

# Clean rules.
clean:
	rm -f $(LIB_TARGET) 
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
