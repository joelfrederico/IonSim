#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "fields.h"
#include "support_func.h"

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
	field(X_PTS, Y_PTS, X_EDGE_MAG, Y_EDGE_MAG),
	x_edge_mag(X_EDGE_MAG),
	y_edge_mag(Y_EDGE_MAG),
	z_edge_mag(Z_EDGE_MAG)
{
}

void FieldTest::custom_init()
{
	int count = 0;
	for (int j=0; j < y_size; j++)
	{
		for (int i=0; i < x_size; i++)
		{
			field.x(i, j) = count;
			field.y(i, j) = -count;
			count++;
		}
	}
}


TEST_F(FieldTest, IsZeroInitially)
{
	bool pass = true;
	double re, im;

	for (long i=0; i < x_size; i++)
	{
		for (long j=0; j < y_size; j++)
		{
			EXPECT_EQ(field.x(i, j), 0) << "At location (" << i << ", " << j <<")";
			EXPECT_EQ(field.y(i, j), 0) << "At location (" << i << ", " << j <<")";
		}
	}
}

TEST_F(FieldTest, AssigmentWorks)
{
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			field.x(i, j) = i*j;
			field.y(i, j) = i*j*2;
		}
	}

	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			EXPECT_EQ(field.x(i, j), i*j)   << "At location (" << i << ", " << j << ")";
			EXPECT_EQ(field.y(i, j), i*j*2) << "At location (" << i << ", " << j << ")";
		}
	}
}

TEST_F(FieldTest, IndicesAreGood)
{
	custom_init();
	for (int i=0; i < x_size*y_size; i++)
	{
		EXPECT_EQ(field.x_data[i], i);
		EXPECT_EQ(field.y_data[i], -i);
	}
}

TEST_F(FieldTest, ArrayIsGood)
{
	custom_init();
	double ** arr;
	int rowCount;
	rowCount = field.x_array_alloc(arr, 0);
	for (long j=0; j < y_size; j++)
	{
		for (long i=0; i < x_size; i++)
		{
			EXPECT_EQ(field.x(i, j), arr[i][j]) << "At location (" << i << ", " << j << ")";
		}
	}
	ionsim::dealloc_2d_array(arr, rowCount);
}

TEST(IonsimTest, AllocThenDealloc)
{
	double ** arr = ionsim::alloc_2d_array(X_PTS, Y_PTS);
	ionsim::dealloc_2d_array(arr, X_PTS);

	arr = ionsim::alloc_2d_array(Y_PTS, X_PTS);
	ionsim::dealloc_2d_array(arr, Y_PTS);
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
