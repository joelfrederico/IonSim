#ifndef __MYGTEST_H_INCLUDED__
#define __MYGTEST_H_INCLUDED__

#include "gtest/gtest.h"
#include "field_data.h"
#include "consts.h"
#include <hdf5.h>

class FieldTest : public ::testing::Test
{
	public:
		FieldTest();

		void custom_init();

		const long x_size;
		const long y_size;
		const long z_size;

		const double x_edge_mag;
		const double y_edge_mag;
		const double z_edge_mag;

		Field_Data field;
};

class HDF5Test : public :: testing::Test
{
	public:
		hid_t file_id;

		HDF5Test();
		~HDF5Test();
};

#endif
