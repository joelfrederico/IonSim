#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "field_data.h"
#include "support_func.h"
#include <gsl/gsl_interp.h>

const long X_PTS = 13;
const long Y_PTS = 15;
const long Z_PTS = 17;

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

void FieldTest::custom_init()
{
	int count = 0;
	for (long k=0; k < z_size; k++)
	{
		for (long j=0; j < y_size; j++)
		{
			for (long i=0; i < x_size; i++)
			{
				field.Ex_ind(i, j, k) = count;
				field.Ey_ind(i, j, k) = -count;
				field.Ez_ind(i, j, k) = count*count;
				count++;
			}
		}
	}
}

TEST_F(FieldTest, DynamicAllocation)
{
	Field_Data *lies;
	lies = new Field_Data(field);
	delete lies;
}

TEST_F(FieldTest, IsZeroInitially)
{
	for (long i=0; i < x_size; i++)
	{
		for (long j=0; j < y_size; j++)
		{
			for (long k=0; k < z_size; k++)
			{
				EXPECT_EQ(field.Ex_ind(i, j, k), 0) << "At location i=" << i << ", j=" << j << ", k=" << k;
				EXPECT_EQ(field.Ey_ind(i, j, k), 0) << "At location i=" << i << ", j=" << j << ", k=" << k;
				EXPECT_EQ(field.Ez_ind(i, j, k), 0) << "At location i=" << i << ", j=" << j << ", k=" << k;
			}
		}
	}
}

TEST_F(FieldTest, AssigmentWorks)
{
	for (long i=0; i < x_size; i++)
	{
		for (long j=0; j < y_size; j++)
		{
			for (long k=0; k < z_size; k++)
			{
				field.Ex_ind(i, j, k) = i*j*k;
				field.Ey_ind(i, j, k) = i*j*k*2;
				field.Ez_ind(i, j, k) = i*j*k*3;
			}
		}
	}

	for (long i=0; i < x_size; i++)
	{
		for (long j=0; j < y_size; j++)
		{
			for (long k=0; k < z_size; k++)
			{
				EXPECT_EQ(field.Ex_ind(i, j, k), i*j*k)   << "At location i=" << i << ", k=" << j << ", k=" << k;
				EXPECT_EQ(field.Ey_ind(i, j, k), i*j*k*2) << "At location i=" << i << ", k=" << j << ", k=" << k;
				EXPECT_EQ(field.Ez_ind(i, j, k), i*j*k*3) << "At location i=" << i << ", k=" << j << ", k=" << k;
			}
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
		EXPECT_EQ(field.z_data[i], i*i);
	}
}

/*
TEST_F(FieldTest, InterpWorks)
{
	custom_init();
	double Ex, Ey, dEx, dEy;

	Ex = field.Ex(13e-6, -5e-6);
	Ey = field.Ey(13e-6, -5e-6);

	dEx = std::abs(84.4333 - Ex);
	dEy = std::abs(-84.4333 - Ey);

	EXPECT_LT(dEx, 4e-5);
	EXPECT_LT(dEy, 4e-5);
}
*/

TEST_F(FieldTest, CopyWorks)
{
	custom_init();
	Field_Data tempfield = field;
	for (int i=0; i < field.x_pts*field.y_pts; i++)
	{
		EXPECT_EQ(tempfield.x_data[i], field.x_data[i]);
		EXPECT_EQ(tempfield.y_data[i], field.y_data[i]);
	}
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
