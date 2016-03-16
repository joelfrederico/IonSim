#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "fields.h"

FieldTest::FieldTest() : x_size(3), y_size(5), z_size(7), field(3, 5, 7)
{
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			for (int k=0; k < z_size; k++)
			{
				field.x(i, j, k) = i*j*k;
				field.y(i, j, k) = i*j*k*2;
			}
		}
	}
}

TEST_F(FieldTest, IsZeroInitially)
{
	long x_size = 3;
	long y_size = 5;
	long z_size = 7;

	bool pass = true;
	double re, im;

	Field field(x_size, y_size, z_size);
	for (long i=0; i < x_size; i++)
	{
		for (long j=0; j < y_size; j++)
		{
			for (long k=0; k < y_size; k++)
			{
				EXPECT_EQ(field.x(i, j, k), 0) << "At location (" << i << ", " << j << ", " << k <<")";
				EXPECT_EQ(field.y(i, j, k), 0) << "At location (" << i << ", " << j << ", " << k <<")";
			}
		}
	}
}

TEST_F(FieldTest, AssigmentWorks)
{
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			for (long k=0; k < y_size; k++)
			{
				EXPECT_EQ(field.x(i, j, k), i*j*k) << "At location (" << i << ", " << j << ", " << k <<")";
				EXPECT_EQ(field.y(i, j, k), i*j*k*2) << "At location (" << i << ", " << j << ", " << k <<")";
			}
		}
	}
}

TEST_F(FieldTest, IndicesAreGood)
{
	Field field(x_size, y_size, z_size);
	int count = 0;
	for (int k=0; k < z_size; k++)
	{
		for (int j=0; j < y_size; j++)
		{
			for (int i=0; i < x_size; i++)
			{
				field.x(i, j, k) = count;
				field.y(i, j, k) = -count;
				count++;
			}
		}
	}

	for (int i=0; i < x_size*y_size*z_size; i++)
	{
		EXPECT_EQ(field.x_data[i], i);
		EXPECT_EQ(field.y_data[i], -i);
	}
}



int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
