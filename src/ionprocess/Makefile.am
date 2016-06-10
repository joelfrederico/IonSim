# AM_CPPFLAGS=$(GSL_CFLAGS) $(HDF5_CPPFLAGS)
AM_CPPFLAGS=$(MYAM_CPPFLAGS)

bin_PROGRAMS = ionprocess

ionprocess_SOURCES=\
		   main.cpp \
		   ../hdf5_classes.cpp

ionprocess_LDADD = $(MYLDADD)
