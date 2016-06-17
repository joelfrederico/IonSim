#include <ionsim.h>
#include "support_func.h"
#include "hdf5test.h"
#include "gtest/gtest.h"

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
