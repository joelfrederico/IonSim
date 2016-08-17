SUBDIRS = ioncheck iontest
AM_CPPFLAGS=$(MYAM_CPPFLAGS) -DPKGDATADIR=\"$(pkgdatadir)\"

bin_PROGRAMS = ionsim iontemp

ionsim_SOURCES = \
		 faddeeva/Faddeeva.cc    \
		 pugixml/src/pugixml.cpp \
		 beam.cpp                \
		 consts.cpp              \
		 ebeam.cpp               \
		 emit.cpp                \
		 fftw_classes.cpp        \
		 field_comm.cpp          \
		 field_data.cpp          \
		 field_interp.cpp        \
		 gsl_classes.cpp         \
		 ions.cpp                \
		 loop_comm.cpp           \
		 main.cpp                \
		 master.cpp              \
		 master_loop.cpp         \
		 mpi_vec.cpp             \
		 parts.cpp               \
		 scalar_data.cpp         \
		 scalar_data_comm.cpp    \
		 simparams.cpp           \
		 slave.cpp               \
		 slave_loops.cpp         \
		 support_func.cpp        \
		 writer_base.cpp         \
		 writer_parallel.cpp     \
		 writer_serial.cpp

ionsim_LDADD = $(MYLDADD)

pkgdata_DATA = ../data/config.xml

iontemp_SOURCES = \
		  pugixml/src/pugixml.cpp \
		  beam.cpp                \
		  emit.cpp                \
		  iontemp.cpp             \
		  simparams.cpp           \
		  support_func.cpp

iontemp_LDADD = $(MYLDADD)
