#ifndef __HDF5_CLASSES_H_INCLUDED__
#define __HDF5_CLASSES_H_INCLUDED__

#include <hdf5.h>
#include <string>
#include <typeinfo>
#include <iostream>

class DataspaceCreate;

class Debug
{
	protected:
		const bool debug;

	public:
		Debug(const bool debug);
};

class GroupAccess : public Debug
{
	protected:
		hid_t _loc_id;
		const std::string _group_str;

	public:
		hid_t group_id;

		GroupAccess(hid_t &loc_id, std::string group_str);
		~GroupAccess();

};

class GroupStepAccess : public GroupAccess
{
	public:
		GroupStepAccess(hid_t &loc_id, unsigned int step);
};

class DatasetAccess : public Debug
{
	private:
		hid_t _loc_id;
		DataspaceCreate *dataspace;

		const std::string _dataset_str;

	public:
		hid_t dataspace_id;
		hid_t dataset_id;

		DatasetAccess(hid_t &loc_id, std::string dataset_str, int rank, hsize_t *count);
		~DatasetAccess();
};

class DataspaceCreate : public Debug
{
	public:
		hid_t dataspace_id;

		DataspaceCreate(int rank, hsize_t *count);
		DataspaceCreate(H5S_class_t type);
		~DataspaceCreate();
};

class PlistCreate : public Debug
{
	private:
		hid_t _cls_id;

	public:
		hid_t plist_id;

		PlistCreate(hid_t cls_id);
		~PlistCreate();
};

class AttributeCreate : public Debug
{
	private:
		hid_t _loc_id;
		hid_t type_id;
		const std::string _attr_name;
	public:
		hid_t attr_id;

		template <class T>
		AttributeCreate(hid_t loc_id, const std::string &attr_name, T attr_value) :
			Debug(true)
		{
			/* hid_t type_id = H5T_NATIVE_DOUBLE; */
			if (typeid(attr_value) == typeid(double)) {
				type_id = H5T_NATIVE_DOUBLE;
			} else if (typeid(attr_value) == typeid(int)) {
				type_id = H5T_NATIVE_INT;
			} else if (typeid(attr_value) == typeid(long)) {
				type_id = H5T_NATIVE_LONG;
			} else if (typeid(attr_value) == typeid(unsigned int)) {
				type_id = H5T_NATIVE_UINT;
			} else {
				printf("Not a valid attribute type.\n");
				return;
			}

			herr_t status;
			const hsize_t mysize = 1;

			// ==================================
			// Create an appropriate dataspace
			// ==================================
			/* dataspace_id = H5Screate(H5S_SCALAR); */
			DataspaceCreate dataspace(H5S_SCALAR);
			status = H5Sset_extent_none(dataspace.dataspace_id);

			// ==================================
			// Create and write attribute
			// ==================================
			attr_id = H5Acreate(loc_id, attr_name.c_str(), type_id, dataspace.dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

			if (attr_id < 0) {
				std::cout << "Attribute not created!" << std::endl;
			} else {
				status = H5Awrite(attr_id, type_id, &attr_value);
			}
		}

		~AttributeCreate();
};


#endif
