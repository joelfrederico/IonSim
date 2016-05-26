#ifndef __HDF5TEST_H_INCLUDED__
#define __HDF5TEST_H_INCLUDED__

#include "gtest/gtest.h"
#include <hdf5.h>

class HDF5Test : public :: testing::Test
{
	public:
		hid_t file_id;

		HDF5Test();
		~HDF5Test();
};

#endif
