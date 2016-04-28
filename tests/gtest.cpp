#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "field_data.h"
#include "support_func.h"
#include <gsl/gsl_interp.h>
#include "hdf5_classes.h"

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

HDF5Test::HDF5Test()
{
	file_id = H5Fcreate("test.h5", H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
}

HDF5Test::~HDF5Test()
{
	herr_t status;

	status = H5Fclose(file_id);
	
	if (status >= 0)
	{
		remove("test.h5");
	}
}

TEST_F(HDF5Test, CreateWorked)
{
	EXPECT_GE(file_id, 0);
}

TEST_F(HDF5Test, GroupCreate)
{
	GroupAccess group(file_id, "testgroup");
	EXPECT_GE(group.group_id, 0);
}

TEST_F(HDF5Test, GroupAccess)
{
	GroupAccess *group1;
	
	group1 = new GroupAccess(file_id, "testgroup");
	delete group1;

	GroupAccess group2(file_id, "testgroup");
	
	EXPECT_GE(group2.group_id, 0);
}

TEST_F(HDF5Test, GroupStepAccess)
{
	H5O_info_t info;

	GroupStepAccess group1(file_id, 0);

	EXPECT_GE(group1.group_id, 0);
}

TEST_F(HDF5Test, DatasetAccess)
{
	int rank = 2;

	hsize_t *count;
	count = new hsize_t[2];
	count[0] = 10;
	count[1] = 6;

	DatasetAccess dataset(file_id, "test", 2, count);

	delete count;

	EXPECT_GE(dataset.dataspace_id, 0);
	EXPECT_GE(dataset.dataset_id, 0);
}

TEST_F(HDF5Test, DataspaceCreate)
{
	int rank = 2;

	hsize_t *count;
	count = new hsize_t[2];
	count[0] = 10;
	count[1] = 6;

	DatasetAccess dataspace(file_id, "test", rank, count);

	delete count;

	EXPECT_GE(dataspace.dataspace_id, 0);
}

TEST_F(HDF5Test, PlistCreate)
{
	PlistCreate plist(H5P_DATASET_XFER);

	EXPECT_GE(plist.plist_id, 0);
}

TEST_F(HDF5Test, AttributeCreateDouble)
{
	herr_t status;
	double a_in = 1.3245;
	double a_out;

	AttributeCreate attr(file_id, "double", a_in);

	status = H5Aread(attr.attr_id, H5T_NATIVE_DOUBLE, &a_out);

	EXPECT_GE(status, 0);
	EXPECT_EQ(a_in, a_out);
}

TEST_F(HDF5Test, AttributeCreateInt)
{
	herr_t status;
	int a_in = 32767;
	int a_out;

	AttributeCreate attr(file_id, "int", a_in);

	status = H5Aread(attr.attr_id, H5T_NATIVE_INT, &a_out);

	EXPECT_GE(status, 0);
	EXPECT_EQ(32767, a_out);
}

TEST_F(HDF5Test, AttributeCreateLong)
{
	herr_t status;
	long a_in = 2147483647;
	long a_out;

	AttributeCreate attr(file_id, "long", a_in);

	status = H5Aread(attr.attr_id, H5T_NATIVE_LONG, &a_out);

	EXPECT_GE(status, 0);
	EXPECT_EQ(2147483647, a_out);
}

TEST_F(HDF5Test, AttributeCreateUnsignedInt)
{
	herr_t status;
	unsigned int a_in = 65535;
	unsigned int a_out;

	AttributeCreate attr(file_id, "long", a_in);

	status = H5Aread(attr.attr_id, H5T_NATIVE_UINT, &a_out);

	EXPECT_GE(status, 0);
	EXPECT_EQ(65535, a_out);
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
