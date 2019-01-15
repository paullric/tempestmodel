# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
# Flags for KNL compilation
CXX=               icpc
F90=               ifort
MPICXX=            mpiicpc
MPIF90=            mpiifort
CXXFLAGS+=         -xMIC-AVX512 -fPIC -g -debug inline-debug-info -qopt-report=3 -I $(SNIPER_ROOT)/include
F90FLAGS+=         -xMIC-AVX512 -fPIC -g -debug inline-debug-info -qopt-report=3
F90_RUNTIME=       -lgfortran
NETCDF=            FALSE
PETSC=             FALSE
LAPACK_INTERFACE=  FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=
LAPACK_LDFLAGS=    -mkl=sequential
# DO NOT DELETE
