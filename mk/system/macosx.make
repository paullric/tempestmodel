# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Mac OS X (Paul's Laptop)

CXX=               g++
F90=               gfortran
MPICXX=            mpic++
MPIF90=            mpif90

F90_RUNTIME=       -L/usr/local/lib -lgfortran

# NetCDF
NETCDF_ROOT=       /opt/local
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib

# PetSc
PETSC_ROOT=        /opt/local/lib/petsc
X11_ROOT=          /opt/X11
PETSC_CXXFLAGS=    -I$(PETSC_ROOT)/include
PETSC_LIBRARIES=   -lpetsc -lX11
PETSC_LDFLAGS=     -L$(PETSC_ROOT)/lib -L$(X11_ROOT)/lib

# LAPACK (Mac OS X Accelerate Framework)
LAPACK_INTERFACE=  FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=  
LAPACK_LDFLAGS=    -framework accelerate

# DO NOT DELETE
