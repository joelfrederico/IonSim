#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "fields.h"

FieldTest::FieldTest() : x_size(10), y_size(10), field(10, 10)
{
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			field.x(i, j) = i*j;
			field.y(i, j) = i*j*2;
		}
	}
}

TEST_F(FieldTest, IsZeroInitially)
{
	long x_size = 10;
	long y_size = 10;
	complex_double test;
	bool pass = true;
	double re, im;

	Field field(x_size, y_size);
	for (long i=0; i < x_size; i++)
	{
		for (long j=0; j < y_size; j++)
		{
			EXPECT_EQ(field.x(i, j), 0) << "At location (" << i << ", " << j << ")";
			EXPECT_EQ(field.y(i, j), 0) << "At location (" << i << ", " << j << ")";
		}
	}
}

TEST_F(FieldTest, AssigmentWorks)
{
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			EXPECT_EQ(field.x(i, j), i*j) << "At location (" << i << ", " << j << ")";
			EXPECT_EQ(field.y(i, j), i*j*2) << "At location (" << i << ", " << j << ")";
		}
	}
}



int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
