# AM_CPPFLAGS=$(GSL_CFLAGS) $(HDF5_CPPFLAGS)
AM_CPPFLAGS=$(MYAM_CPPFLAGS)

bin_PROGRAMS = iontest

iontest_SOURCES=\
		../pugixml/src/pugixml.cpp \
		 main.cpp \
		 ../ions.cpp \
		 ../field_interp.cpp \
		 ../field_data.cpp \
		 ../parts.cpp \
		 ../beam.cpp \
		 ../emit.cpp \
		 ../support_func.cpp \
		 ../simparams.cpp

iontest_LDADD = $(MYLDADD)

iontest_CPPFLAGS = -I$(top_srcdir)/src $(MYAM_CPPFLAGS)
