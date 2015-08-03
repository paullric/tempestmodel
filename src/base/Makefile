# TempestBase directory
TEMPESTBASEDIR= ../..

# Library name
LIBNAME= libtempestbase

# Use BLAS
USEBLAS= true

# Load system-specific defaults
include $(TEMPESTBASEDIR)/mk/Make.defs

# Files in this directory
FILES= Preferences.cpp \
	DataContainer.cpp \
	FunctionTimer.cpp \
	MathHelper.cpp \
	Exception.cpp \
	Announce.cpp \
	LinearAlgebra.cpp \
	LegendrePolynomial.cpp \
	PolynomialInterp.cpp \
	MemoryTools.cpp \
	GaussQuadrature.cpp \
	GaussLobattoQuadrature.cpp \
	TimeObj.cpp

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
