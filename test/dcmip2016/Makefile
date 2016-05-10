# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTBASEDIR= ../..

# Load Makefile framework. 
include $(TEMPESTBASEDIR)/mk/framework.make

FILES= TropicalCycloneTest.cpp \
       BaroclinicWaveUMJSTest.cpp \
       SupercellTest.cpp

SHARED_FILES= KesslerPhysics.cpp \
              DCMIPPhysics.cpp

INTERFACE_FILES= interface/tropical_cyclone_test.f90 \
                 interface/baroclinic_wave_test.f90 \
                 interface/supercell_test.f90 \
                 interface/kessler.f90 \
                 interface/simple_physics_v6.f90

SHARED_OBJS= $(SHARED_FILES:%.cpp=$(BUILDDIR)/%.o) 
INTERFACE_OBJS= $(INTERFACE_FILES:interface/%.f90=$(BUILDDIR)/interface/%.o) 

EXEC_TARGETS= $(FILES:%.cpp=%)
CLEAN_TARGETS= $(addsuffix .clean,$(EXEC_TARGETS))

.PHONY: all clean

# Build rules. 
all: $(EXEC_TARGETS)

$(EXEC_TARGETS): %: $(BUILDDIR)/%.o $(SHARED_OBJS) $(INTERFACE_OBJS) $(TEMPESTLIBS)
	$(CXX) $(LDFLAGS) -o $@ $(BUILDDIR)/$*.o $(SHARED_OBJS) $(INTERFACE_OBJS) $(LIBRARIES) $(F90_RUNTIME)

# Go up to the root directory and build the libraries if needed.
$(TEMPESTBASEDIR)/%.a:
	cd $(TEMPESTBASEDIR); $(MAKE) $*.a

# Clean rules.
clean: $(CLEAN_TARGETS)
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)
	rm -f *.mod

$(CLEAN_TARGETS): %.clean:
	rm -f $*

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)
-include $(SHARED_FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
