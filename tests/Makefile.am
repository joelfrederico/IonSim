# Bulid rules for samples and tests. Automake's naming for some of
# these variables isn't terribly obvious, so this is a brief
# reference:
#
# TESTS -- Programs run automatically by "make check"
# check_PROGRAMS -- Programs built by "make check" but not necessarily run
AM_CPPFLAGS=$(MYAM_CPPFLAGS)

check_LTLIBRARIES = libgtest.la
libgtest_la_SOURCES = ../googletest/googletest/src/gtest-all.cc
libgtest_la_CPPFLAGS =-I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest $(MYAM_CPPFLAGS)
libgtest_la_LDFLAGS =-pthread

check_PROGRAMS = gtest mpi_gtest

gtest_SOURCES = gtest.cpp \
		fieldtest.cpp \
		hdf5test.cpp \
		ionstest.cpp \
		simparamstest.cpp \
		support.cpp \
		../src/field_data.cpp \
		../src/field_comm.cpp \
		../src/hdf5_classes.cpp \
		../src/ions.cpp \
		../src/loop_comm.cpp \
		../src/support_func.cpp \
		../src/simparams.cpp \
		../src/emit.cpp \
		../src/beam.cpp

gtest_LDADD = libgtest.la $(MYLDADD)

gtest_LDFLAGS = -pthread -L/opt/local/lib

gtest_CPPFLAGS =-I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/src $(MYAM_CPPFLAGS)

mpi_gtest_SOURCES = mpi_gtest.cpp

mpi_gtest_LDADD = libgtest.la $(MYLDADD)
mpi_gtest_LDFLAGS = -pthread -L/opt/local/lib
mpi_gtest_CPPFLAGS =-I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/src $(MYAM_CPPFLAGS)
 
TESTS = gtest
