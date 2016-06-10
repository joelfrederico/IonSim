#include "hdf5_classes.h"
#include <typeinfo>

// ==================================
// Debug
// ==================================
Debug::Debug(const bool debug_flag) :
	debug(debug_flag)
{
}

herr_t Debug::close(herr_t (*f)(hid_t), hid_t attr_id, std::string name)
{
	herr_t status;

	status = f(attr_id);
	if (debug)
	{
		if (status < 0)
		{
			std::cout << "Couldn't close: " << name << std::endl;
		} else {
			std::cout << "Successfully closed: " << name << std::endl;
		}
	}

	return status;
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
		if ((status < 0) || (group_id < 0))
		{
			std::cout << "Warning: could not create " << _group_str << std::endl;
		}
	}
}

GroupAccess::~GroupAccess()
{
	/* herr_t status; */

	/* status = H5Gclose(group_id); */

	close(&H5Gclose, group_id, _group_str);
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
// DatasetOpen
// ==================================
DatasetOpen::DatasetOpen(hid_t &loc_id, std::string dataset_str) :
	Debug(false),
	_dataset_str(dataset_str)
{
	_loc_id = loc_id;
	
	// ==================================
	// Open dataset
	// ==================================
	dataset_id = H5Dopen(_loc_id, dataset_str.c_str(), H5P_DEFAULT);
}

DatasetOpen::~DatasetOpen()
{
	close(&H5Dclose, dataset_id, _dataset_str);
}

// ==================================
// DatasetAccess
// ==================================
DatasetAccess::DatasetAccess(hid_t &loc_id, std::string dataset_str, int rank, hsize_t *count) :
	Debug(false),
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
	close(&H5Dclose, dataset_id, _dataset_str);

	delete dataspace;
}

// ==================================
// DataspaceCreate
// ==================================
DataspaceCreate::DataspaceCreate(int rank, hsize_t *count) :
	Debug(false)
{
	dataspace_id = H5Screate_simple(rank, count, NULL);
}

DataspaceCreate::DataspaceCreate(H5S_class_t type) :
	Debug(false)
{
	dataspace_id = H5Screate(type);
}

DataspaceCreate::~DataspaceCreate()
{
	close(&H5Sclose, dataspace_id, "dataspace");
}

// ==================================
// PlistCreate
// ==================================
PlistCreate::PlistCreate(hid_t cls_id) :
	Debug(false)
{
	_cls_id = cls_id;
	plist_id = H5Pcreate(cls_id);
}

PlistCreate::~PlistCreate()
{
	H5Pclose(plist_id);
}

// ==================================
// AttributeOpen
// ==================================
AttributeOpen::AttributeOpen(hid_t loc_id, std::string attr_name) :
	Debug(false),
	_loc_id(loc_id)
{
	attr_id = H5Aopen(loc_id, attr_name.c_str(), H5P_DEFAULT);
}

AttributeOpen::~AttributeOpen()
{
	close(&H5Aclose, attr_id, _attr_name);
}

int AttributeOpen::read()
{
	int out;
	herr_t status;

	status = H5Aread(attr_id, H5T_NATIVE_INT, &out);

	return out;
}

// ==================================
// AttributeCreate
// ==================================
AttributeCreate::~AttributeCreate()
{
	close(&H5Aclose, attr_id, _attr_name);

	delete _dataspace;
}

// ==================================
// FileOpen
// ==================================
hid_t _fileopen(std::string filename, unsigned flags)
{
	return H5Fopen(filename.c_str(), flags, H5P_DEFAULT);
}

FileOpen::FileOpen(std::string filename) :
	Debug(false)
{
	_filename = filename;
	file_id = _fileopen(filename, H5F_ACC_RDONLY);
}

FileOpen::FileOpen(std::string filename, unsigned flags) :
	Debug(false)
{
	_filename = filename;
	file_id = _fileopen(filename, flags);
}

FileOpen::~FileOpen()
{
	close(&H5Fclose, file_id, _filename);
}

// ==================================
// FileCreate
// ==================================
hid_t _filecreate(std::string filename, unsigned flags)
{
	return H5Fcreate(filename.c_str(), flags, H5P_DEFAULT, H5P_DEFAULT);
}

FileCreate::FileCreate(std::string filename) :
	Debug(false)
{
	_filename = filename;
	file_id = _filecreate(filename, H5F_ACC_TRUNC);
}

FileCreate::FileCreate(std::string filename, unsigned flags) :
	Debug(false)
{
	_filename = filename;
	file_id = _filecreate(filename, flags);
}

FileCreate::~FileCreate()
{
	close(&H5Fclose, file_id, _filename);
}
