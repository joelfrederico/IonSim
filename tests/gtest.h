#ifndef __MYGTEST_H_INCLUDED__
#define __MYGTEST_H_INCLUDED__

#include "gtest/gtest.h"
#include "fields.h"
#include "consts.h"

class FieldTest : public ::testing::Test
{
	public:
		FieldTest();

		long x_size;
		long y_size;
		long z_size;

		Field field;
};

#endif
