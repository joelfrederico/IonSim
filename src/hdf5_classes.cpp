#include "hdf5_classes.h"
#include <typeinfo>

// ==================================
// Debug
// ==================================
Debug::Debug(const bool debug_flag) :
	debug(debug_flag)
{
}

// ==================================
// GroupAccess
// ==================================
GroupAccess::GroupAccess(hid_t &loc_id, std::string group_str) :
	Debug(false),
	_group_str(group_str)
{
	_loc_id = loc_id;

	herr_t status;
	H5G_info_t objinfo;

	// ==================================
	// Access or create a new group
	// ==================================
	H5Eset_auto(NULL, NULL, NULL);
	group_id = H5Gopen(loc_id, _group_str.c_str(), H5P_DEFAULT);
	if (group_id < 0)
	{
		group_id = H5Gcreate(loc_id, _group_str.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gget_info_by_name(loc_id, _group_str.c_str(), &objinfo, H5P_DEFAULT);
		if ((status < 0) || (group_id < 0)) printf("WARNING\n");
	}
}

GroupAccess::~GroupAccess()
{
	herr_t status;

	status = H5Gclose(group_id);

	if (debug)
	{
		if (status < 0)
		{
			std::cout << "Failed to close group: " << _group_str << std::endl;
		} else {
			std::cout << "Successfully closed: " << _group_str << std::endl;
		}
	}
}

// ==================================
// GroupStepAccess
// ==================================
std::string _getstep(unsigned int step)
{
	std::string out;
	char buf[10];
	sprintf(buf, "Step_%04u", step);
	out = buf;
	return out;
}

GroupStepAccess::GroupStepAccess(hid_t &loc_id, unsigned int step) :
	GroupAccess(loc_id, _getstep(step))
{
}

// ==================================
// DatasetAccess
// ==================================
DatasetAccess::DatasetAccess(hid_t &loc_id, std::string dataset_str, int rank, hsize_t *count) :
	Debug(true),
	_dataset_str(dataset_str)
{
	_loc_id = loc_id;
	
	dataspace = new DataspaceCreate(rank, count);
	dataspace_id = (*dataspace).dataspace_id;

	// ==================================
	// Create dataset
	// ==================================
	dataset_id = H5Dcreate(_loc_id, dataset_str.c_str(), H5T_NATIVE_DOUBLE, (*dataspace).dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

DatasetAccess::~DatasetAccess()
{
	H5Dclose(dataset_id);
	delete dataspace;
}

// ==================================
// DataspaceCreate
// ==================================
DataspaceCreate::DataspaceCreate(int rank, hsize_t *count) :
	Debug(true)
{
	dataspace_id = H5Screate_simple(rank, count, NULL);
}

DataspaceCreate::DataspaceCreate(H5S_class_t type) :
	Debug(true)
{
	dataspace_id = H5Screate(type);
}

DataspaceCreate::~DataspaceCreate()
{
	H5Sclose(dataspace_id);
}

// ==================================
// PlistCreate
// ==================================
PlistCreate::PlistCreate(hid_t cls_id) :
	Debug(true)
{
	_cls_id = cls_id;
	plist_id = H5Pcreate(cls_id);
}

PlistCreate::~PlistCreate()
{
	H5Pclose(plist_id);
}

// ==================================
// AttributeCreate
// ==================================
AttributeCreate::~AttributeCreate()
{
	H5Aclose(attr_id);
	delete _dataspace;
}
