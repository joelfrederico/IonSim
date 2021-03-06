#ifndef __FIELDTEST_H_INCLUDED__
#define __FIELDTEST_H_INCLUDED__

#include "gtest/gtest.h"
#include "field_data.h"

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

#endif
