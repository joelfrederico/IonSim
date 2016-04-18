SUBDIRS = ioncheck
AM_CPPFLAGS=$(MYAM_CPPFLAGS)

bin_PROGRAMS=ionsim

ionsim_SOURCES=\
		 beam.cpp             \
		 consts.cpp           \
		 ebeam.cpp            \
		 emit.cpp             \
		 faddeeva/Faddeeva.cc \
		 field_data.cpp       \
		 field_interp.cpp     \
		 field_comm.cpp       \
		 hdf5_classes.cpp     \
		 ions.cpp             \
		 loop_comm.cpp        \
		 main.cpp             \
		 master.cpp           \
		 parts.cpp            \
		 simparams.cpp        \
		 slave.cpp            \
		 support_func.cpp     \
		 writer_base.cpp      \
		 writer_serial.cpp    \
		 writer_parallel.cpp

ionsim_LDADD = $(MYLDADD)
