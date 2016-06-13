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

		herr_t close(herr_t (*f)(hid_t), hid_t attr_id, std::string name);
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

class DatasetOpen : public Debug
{
	private:
		hid_t _loc_id;
		/* DataspaceCreate *dataspace; */

	public:
		const std::string _dataset_str;

		/* hid_t dataspace_id; */
		hid_t dataset_id;
		hid_t dataspace_id;

		DatasetOpen(hid_t &loc_id, std::string dataset_str);
		~DatasetOpen();
};

class DatasetAccess : public Debug
{
	private:
		hid_t _loc_id;
		DataspaceCreate *dataspace;

	public:
		const std::string _dataset_str;

		hid_t dataspace_id;
		hid_t dataset_id;

		DatasetAccess(hid_t &loc_id, std::string dataset_str, int rank, hsize_t *dims);
		~DatasetAccess();
};

class DataspaceCreate : public Debug
{
	public:
		hid_t dataspace_id;

		DataspaceCreate(int rank, hsize_t *dims);
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

class AttributeOpen : public Debug
{
	private:
		const std::string _attr_name;
		const hid_t _loc_id;

	public:
		hid_t attr_id;

		AttributeOpen(hid_t loc_id, std::string attr_name);
		~AttributeOpen();

		int read();
};

class AttributeCreate : public Debug
{
	private:
		hid_t _loc_id;
		hid_t type_id;
		const std::string _attr_name;
		DataspaceCreate * _dataspace;

	public:
		hid_t attr_id;

		template <class T>
		AttributeCreate(hid_t loc_id, const std::string &attr_name, T attr_value) :
			Debug(false),
			_attr_name(attr_name)
		{
			herr_t status;

			// ==================================
			// Figure out the HDF5 type to write
			// ==================================
			if (typeid(attr_value) == typeid(double)) {
				type_id = H5T_NATIVE_DOUBLE;
			} else if (typeid(attr_value) == typeid(int)) {
				type_id = H5T_NATIVE_INT;
			} else if (typeid(attr_value) == typeid(long)) {
				type_id = H5T_NATIVE_LONG;
			} else if (typeid(attr_value) == typeid(long long)) {
				type_id = H5T_NATIVE_LLONG;
			} else if (typeid(attr_value) == typeid(unsigned int)) {
				type_id = H5T_NATIVE_UINT;
			} else if (typeid(attr_value) == typeid(const char *)) {
				type_id = H5Tcopy(H5T_C_S1);
				status = H5Tset_size(type_id, H5T_VARIABLE);
			} else {
				printf("Not a valid attribute type.\n");
				return;
			}

			// ==================================
			// Create an appropriate dataspace
			// ==================================
			_dataspace = new DataspaceCreate(H5S_SCALAR);
			status = H5Sset_extent_none(_dataspace->dataspace_id);

			// ==================================
			// Create and write attribute
			// ==================================
			attr_id = H5Acreate(loc_id, attr_name.c_str(), type_id, _dataspace->dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

			if (attr_id < 0) {
				std::cout << "Attribute not created!" << std::endl;
			} else {
				status = H5Awrite(attr_id, type_id, &attr_value);
			}

			// ==================================
			// Close resources
			// ==================================
			if (typeid(attr_value) == typeid(const char *)) {
				H5Tclose(type_id);
			}
		}

		~AttributeCreate();
};

class FileOpen : public Debug
{
	private:
		std::string _filename;

	public:
		hid_t file_id;

		FileOpen(std::string filename);
		FileOpen(std::string filename, unsigned flags);
		~FileOpen();
};

class FileCreate : public Debug
{
	private:
		std::string _filename;

	public:
		hid_t file_id;

		FileCreate(std::string filename);
		FileCreate(std::string filename, unsigned flags);
		~FileCreate();
};

#endif
