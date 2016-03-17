#ifndef __MYGTEST_H_INCLUDED__
#define __MYGTEST_H_INCLUDED__

#include "gtest/gtest.h"
#include "fields.h"
#include "consts.h"

class FieldTest : public ::testing::Test
{
	public:
		FieldTest();

		const long x_size;
		const long y_size;
		const long z_size;

		const double x_edge_mag;
		const double y_edge_mag;
		const double z_edge_mag;

		Field field;
};

#endif
