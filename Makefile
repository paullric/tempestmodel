# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

LIB_TARGETS= lib/libtempestbase.a lib/libhardcoreatm.a 
CLEAN_TARGETS= src/base.clean src/atm.clean test.clean

.PHONY: all test clean $(CLEAN_TARGETS)

# Build rules.
all: $(LIB_TARGETS) test

lib/libtempestbase.a:
	cd src/base; $(MAKE)

lib/libhardcoreatm.a:
	cd src/atm; $(MAKE)

test: $(LIB_TARGETS) 
	cd test; $(MAKE)

# Clean rules.
clean: $(CLEAN_TARGETS)

$(CLEAN_TARGETS): %.clean:
	cd $*; $(MAKE) clean

# DO NOT DELETE
