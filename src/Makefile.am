SUBDIRS = ioncheck iontest
AM_CPPFLAGS=$(MYAM_CPPFLAGS)

bin_PROGRAMS = ionsim iontemp

ionsim_SOURCES = \
		 faddeeva/Faddeeva.cc    \
		 pugixml/src/pugixml.cpp \
		 beam.cpp                \
		 consts.cpp              \
		 ebeam.cpp               \
		 emit.cpp                \
		 field_comm.cpp          \
		 field_data.cpp          \
		 field_interp.cpp        \
		 hdf5_classes.cpp        \
		 ions.cpp                \
		 loop_comm.cpp           \
		 main.cpp                \
		 master.cpp              \
		 parts.cpp               \
		 simparams.cpp           \
		 slave.cpp               \
		 support_func.cpp        \
		 writer_base.cpp         \
		 writer_parallel.cpp     \
		 writer_serial.cpp

ionsim_LDADD = $(MYLDADD)

iontemp_SOURCES = \
		  pugixml/src/pugixml.cpp \
		  beam.cpp                \
		  emit.cpp                \
		  iontemp.cpp             \
		  simparams.cpp           \
		  support_func.cpp

iontemp_LDADD = $(MYLDADD)
