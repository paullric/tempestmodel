# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTBASEDIR= ../..

# Load Makefile framework. 
include $(TEMPESTBASEDIR)/mk/framework.make

FILES= BaroclinicWaveJWTest.cpp \
       BaroclinicWaveUMJSTest.cpp \
       MountainWaveSphereTest.cpp \
       ScharMountainSphereTest.cpp \
       InertiaGravityWaveTest.cpp \
       BaldaufGravityWaveTest.cpp \
       MountainRossby3DTest.cpp \
       StationaryMountainFlowTest.cpp \
       HeldSuarezTest.cpp

EXEC_TARGETS= $(FILES:%.cpp=%)
CLEAN_TARGETS= $(addsuffix .clean,$(EXEC_TARGETS))

.PHONY: all clean

# Build rules. 
all: $(EXEC_TARGETS)

$(EXEC_TARGETS): %: $(BUILDDIR)/%.o $(TEMPESTLIBS)
	$(CXX) $(LDFLAGS) -o $@ $(BUILDDIR)/$*.o $(LIBRARIES) 

# Go up to the root directory and build the libraries if needed.
$(TEMPESTBASEDIR)/%.a:
	cd $(TEMPESTBASEDIR); $(MAKE) $*.a

# Clean rules.
clean: $(CLEAN_TARGETS)
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

$(CLEAN_TARGETS): %.clean:
	rm -f $*

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
