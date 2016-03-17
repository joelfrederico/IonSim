#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "fields.h"

const long X_PTS = 3;
const long Y_PTS = 5;
const long Z_PTS = 7;

const double X_EDGE_MAG = 30e-6;
const double Y_EDGE_MAG = 30e-6;
const double Z_EDGE_MAG = 30e-6;

FieldTest::FieldTest() : 
	x_size(X_PTS),
	y_size(Y_PTS),
	z_size(Z_PTS),
	field(X_PTS, Y_PTS, Z_PTS, X_EDGE_MAG, Y_EDGE_MAG, Z_EDGE_MAG),
	x_edge_mag(X_EDGE_MAG),
	y_edge_mag(Y_EDGE_MAG),
	z_edge_mag(Z_EDGE_MAG)
{
}

TEST_F(FieldTest, IsZeroInitially)
{
	bool pass = true;
	double re, im;

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
			for (int k=0; k < z_size; k++)
			{
				field.x(i, j, k) = i*j*k;
				field.y(i, j, k) = i*j*k*2;
			}
		}
	}

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
