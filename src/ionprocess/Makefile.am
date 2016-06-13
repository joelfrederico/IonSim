# AM_CPPFLAGS=$(GSL_CFLAGS) $(HDF5_CPPFLAGS)
AM_CPPFLAGS=$(MYAM_CPPFLAGS)

bin_PROGRAMS = ionprocess

lib_LTLIBRARIES = libionprocess.la

include_HEADERS = ionprocess.h

libionprocess_la_SOURCES=\
			 ionprocess.cpp \
			 ../hdf5_classes.cpp

ionprocess_SOURCES=\
		   main.cpp \
		   ../hdf5_classes.cpp

ionprocess_LDADD = libionprocess.la $(MYLDADD)

ionprocess_CPPFLAGS = $(AM_CPPFLAGS)
