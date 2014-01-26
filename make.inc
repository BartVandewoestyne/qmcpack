# To compile the code, run something like:
#
#       $ make FC=g95 MODE=optimized
#
# Possible values for FC are:
#
#       g95, gfortran, F, nagfor, ifc, ifort
#
# Possible values for MODE are:
#
#       debug, optimized, profile
#
# The default compiler and mode are specified here.
#
FC = gfortran
MODE = debug


# The Open Source g95
# Notes:
#   * g95 had no clean gdb support on 15/07/2005.
#   * -O4 is invalid for g95 on 03/10/2005.
#   * -O2 might often be faster than -O3
#   * performance wise, the following options are also sometimes advised
#     for g95:
#       -march=pentium4 -funroll-loops -mfpmath=sse -malign-double
#
FFLAGS_g95_debug = -g -pedantic -Wall -fbounds-check -ftrace=full -std=F
FFLAGS_g95_optimized = -O3 -ffast-math -s
FFLAGS_g95_profile = -pg


# The Open Source gfortran.
# Notes:
#   * gfortran also supports OpenMP by the -fopenmp option.
#   * performancewise, the following options might also be useful:
#       -march=pentium4 -funroll-loops -mfpmath=sse -malign-double
#
FFLAGS_gfortran_debug = -g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none -std=f95
FFLAGS_gfortran_optimized = -O3 -ffast-math
FFLAGS_gfortran_profile = -pg


# The free F-compiler
# 
FFLAGS_F_debug = -g -C=all -gline -info
FFLAGS_F_optimized = -O4
FFLAGS_F_profile = -pg


# The commercial NagWare nagfor compiler
# Notes:
#   * use `upsf95' to debug.
#
FFLAGS_nagfor_debug = -g -C=all -C=undefined -gline -u -info
FFLAGS_nagfor_optimized = -O3 -Oassumed
FFLAGS_nagfor_profile = -pg


# The commercial Intel(R) Fortran Compiler for 32-bit applications, Version 7.0
# Notes:
#   * the -Tf flag must be at the end of FFLAGS
#   * use `idb' to debug
#   * use `prof` (for mon.out) or `gprof` (for gmon.out) for profiling
#
FFLAGS_ifc_debug = -g -e95 -implicitnone -FR -C -Tf
FFLAGS_ifc_optimized = -O2 -FR -Tf
FFLAGS_ifc_profile = -qp


# The commercial Intel(R) Fortran Compiler for 32-bit applications, Version 9.0
# Notes:
#   * the -Tf flag must be at the end of FFLAGS
#   * profile with `gprof'
#
FFLAGS_ifort_debug = -g -stand f95 -implicitnone -check -warn -free -Tf
FFLAGS_ifort_optimized = -O3 -free -Tf
FFLAGS_ifort_profile = -qp

# Sun's Forte compiler
# f95: Forte Developer 7 Fortran 95 7.0 Patch 111714-09 2003/10/15
#FFLAGS_f95_debug = -C -g -u -ansi -free -M../src
#FFLAGS_f95_optimized = -O3 -free -M../src
#FFLAGS_f95_profile = -pg -M../src


FFLAGS = $(FFLAGS_$(FC)_$(MODE))

LOADER = $(FC)

#
# The archiver and the flag(s) to use when building archive (library)
# The 'c' flag tells ar to create the library if it doesn't already exist.
# The 'r' flag tells it to replace older object files in the library with the
# new object files.
#
ARCH      = ar
ARCHFLAGS = cr

# The command used to create or update the index for the archive.
RANLIB = ranlib

#
# The location of the libraries to which you will link.
#
QMCPACKLIB = qmcpack

#
# The directory to install the library
#
INSTALLDIR = $(HOME)/lib
