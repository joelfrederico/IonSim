#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([IonMotion], [0.0.1], [joelfred@slac.stanford.edu])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AX_PROG_CXX_MPI

# Checks for libraries.
PKG_CHECK_MODULES([GSL], [gsl >= 1.16])
AX_LIB_HDF5([parallel])
if test "$with_hdf5" = "yes"; then
	HDF5_LIBS="$HDF5_LIBS -lhdf5_cpp"
fi



# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
