# Library base directory
LIBBASEDIR= ..

# Library name
LIBNAME= libtempestbase

# Use BLAS
USEBLAS= true

# Load system-specific defaults
USEBLAS= true
include $(LIBBASEDIR)/mk/Make.defs

# Files in this directory
FILES= Preferences.cpp \
	MatlabOutput.cpp \
	MatlabInput.cpp \
	FunctionTimer.cpp \
	MathHelper.cpp \
	Announce.cpp \
	GaussianQuadrature.cpp \
	LinearAlgebra.cpp \
	LegendrePolynomial.cpp \
	PolynomialInterp.cpp

##
## Build instructions
##
all: $(FILES:%.cpp=$(BUILDDIR)/%.o)
	cp *.h $(LIBBASEDIR)/include
	ar -cqs $(LIBBASEDIR)/lib/$(LIBNAME).a build/*.o

##
## Clean
##
clean:
	rm -rf build
	rm -rf depend
	rm -rf lib

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
