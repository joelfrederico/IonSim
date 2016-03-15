noinst_PROGRAMS = gtest demo

gtest_SOURCES = gtest.cpp ../src/fields.cpp

gtest_LDADD = libgtest.la

gtest_LDFLAGS = -pthread

gtest_CPPFLAGS = -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/src

# Bulid rules for samples and tests. Automake's naming for some of
# these variables isn't terribly obvious, so this is a brief
# reference:
#
# TESTS -- Programs run automatically by "make check"
# check_PROGRAMS -- Programs built by "make check" but not necessarily run

check_PROGRAMS = gtest

TESTS = $(check_PROGRAMS)

demo_SOURCES = demo.cpp ../src/fields.cpp

demo_CPPFLAGS = -I$(top_srcdir)/src

noinst_LTLIBRARIES = libgtest.la

libgtest_la_SOURCES = ../googletest/googletest/src/gtest-all.cc
libgtest_la_CPPFLAGS = -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread
