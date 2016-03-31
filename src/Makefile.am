# AM_CPPFLAGS=$(GSL_CFLAGS) $(HDF5_CPPFLAGS)
AM_CPPFLAGS=$(MYAM_CPPFLAGS)

bin_PROGRAMS=ionsim

ionsim_SOURCES=\
		 beam.cpp \
		 consts.cpp \
		 ebeam.cpp \
		 emit.cpp \
		 faddeeva/Faddeeva.cc \
		 field_data.cpp \
		 field_comm.cpp \
		 ions.cpp \
		 main.cpp \
		 master.cpp \
		 parts.cpp \
		 simparams.cpp \
		 slave.cpp \
		 support_func.cpp \
		 writer.cpp

ionsim_LDADD = $(MYLDADD)
# ionsim_LDADD = $(GSL_LIBS) $(HDF5_LIBS)
# ionsim_CPPFLAGS = $(GSL_CFLAGS) $(HDF5_CPPFLAGS)
