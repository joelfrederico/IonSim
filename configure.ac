#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([IonMotion], [0.0.1], [joelfred@slac.stanford.edu])
AM_INIT_AUTOMAKE([-Wall -Werror -Woverride subdir-objects])
AC_CONFIG_SRCDIR([src/main.cpp])
# AC_CONFIG_MACRO_DIRS([googletest/googletest/m4])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AX_PROG_CXX_MPI
AM_PROG_AR
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])

# Checks for libraries.
PKG_CHECK_MODULES([GSL], [gsl >= 1.15])
AX_LIB_HDF5([parallel])
if test "$with_hdf5" = "yes"; then
	HDF5_LIBS="$HDF5_LIBS -lhdf5_cpp"
fi

MYLDADD="$GSL_LIBS $HDF5_LIBS"
AC_SUBST(MYLDADD,[$MYLDADD])

MYAM_CPPFLAGS="$GSL_CFLAGS $HDF5_CPPFLAGS"
AC_SUBST(MYAM_CPPFLAGS,[$MYAM_CPPFLAGS])

# AC_LANG_PUSH([C++])
# AC_CHECK_HEADER([gtest/gtest.h], [AC_DEFINE([HAVE_GTEST_H], 1)])
# AC_LANG_POP([C++])
# GTEST_LIB_CHECK([1.7.0])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


AC_CONFIG_FILES([Makefile
                 src/Makefile
		 tests/Makefile
		 ])
AC_OUTPUT
