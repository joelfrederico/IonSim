#include "writer_base.h"
#include "field_data.h"
#include "support_func.h"

WriterBase::WriterBase(const std::string &filename)
{
	_filename = filename;
}

WriterBase::~WriterBase()
{
	H5Fclose(file_id);
}


hid_t WriterBase::group_access(hid_t &loc_id, const std::string &group)
{
	return group_access(loc_id, group.c_str());
}

hid_t WriterBase::group_access(hid_t &loc_id, const char *group)
{
	herr_t status;
	hid_t group_id;
	H5G_info_t objinfo;

	// ==================================
	// Access or create a new group
	// ==================================
	H5Eset_auto(NULL, NULL, NULL);
	group_id = H5Gopen(loc_id, group, H5P_DEFAULT);
	if (group_id < 0)
	{
		group_id = H5Gcreate(loc_id, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gget_info_by_name(loc_id, group, &objinfo, H5P_DEFAULT);
	}
	return group_id;
}

hid_t WriterBase::group_step_access(hid_t &loc_id, long step)
{
	char buf[10];
	sprintf(buf, "Step_%04ld", step);
	return group_access(loc_id, buf);
}

hid_t WriterBase::dataset_create(hid_t &group_id, hid_t &dataspace_id, int rank, hsize_t *count, const std::string &dataset)
{
	hid_t dataset_id;

	dataspace_id = H5Screate_simple(rank, count, NULL);

	// ==================================
	// Create dataset
	// ==================================
	dataset_id = H5Dcreate(group_id, dataset.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	return dataset_id;
}
