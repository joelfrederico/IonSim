#include "fieldtest.h"

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

				field.Bx_ind(i, j, k) = count;
				field.By_ind(i, j, k) = -count;
				field.Bz_ind(i, j, k) = count*count;
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

				EXPECT_EQ(field.Bx_ind(i, j, k), 0) << "At location i=" << i << ", j=" << j << ", k=" << k;
				EXPECT_EQ(field.By_ind(i, j, k), 0) << "At location i=" << i << ", j=" << j << ", k=" << k;
				EXPECT_EQ(field.Bz_ind(i, j, k), 0) << "At location i=" << i << ", j=" << j << ", k=" << k;
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

				field.Bx_ind(i, j, k) = i*j*k;
				field.By_ind(i, j, k) = i*j*k*2;
				field.Bz_ind(i, j, k) = i*j*k*3;
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

				EXPECT_EQ(field.Bx_ind(i, j, k), i*j*k)   << "At location i=" << i << ", k=" << j << ", k=" << k;
				EXPECT_EQ(field.By_ind(i, j, k), i*j*k*2) << "At location i=" << i << ", k=" << j << ", k=" << k;
				EXPECT_EQ(field.Bz_ind(i, j, k), i*j*k*3) << "At location i=" << i << ", k=" << j << ", k=" << k;
			}
		}
	}
}

TEST_F(FieldTest, IndicesAreGood)
{
	custom_init();
	for (int i=0; i < x_size*y_size; i++)
	{
		EXPECT_EQ(field.Ex_data[i], i);
		EXPECT_EQ(field.Ey_data[i], -i);
		EXPECT_EQ(field.Ez_data[i], i*i);

		EXPECT_EQ(field.Bx_data[i], i);
		EXPECT_EQ(field.By_data[i], -i);
		EXPECT_EQ(field.Bz_data[i], i*i);
	}
}

TEST_F(FieldTest, CopyWorks)
{
	custom_init();
	Field_Data tempfield = field;
	for (int i=0; i < field.x_pts*field.y_pts; i++)
	{
		EXPECT_EQ(tempfield.Ex_data[i], field.Ex_data[i]);
		EXPECT_EQ(tempfield.Ey_data[i], field.Ey_data[i]);

		EXPECT_EQ(tempfield.Bx_data[i], field.Bx_data[i]);
		EXPECT_EQ(tempfield.By_data[i], field.By_data[i]);
	}
}
