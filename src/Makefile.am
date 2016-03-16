bin_PROGRAMS = ionsim

ionsim_SOURCES = \
		 baseclass.cpp \
		 beam.cpp \
		 consts.cpp \
		 ebeam.cpp \
		 emit.cpp \
		 faddeeva/Faddeeva.cc \
		 fields.cpp \
		 ions.cpp \
		 main.cpp \
		 master.cpp \
		 simparams.cpp \
		 slave.cpp \
		 support_func.cpp

ionsim_LDADD = $(GSL_LIBS) $(HDF5_LIBS)
ionsim_CPPFLAGS = $(GSL_CFLAGS) $(HDF5_CPPFLAGS)
