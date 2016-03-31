#ifndef __WRITER_H_INCLUDED__
#define __WRITER_H_INCLUDED__

#include <string>
#include <hdf5.h>
#include "field_data.h"
#include "simparams.h"

class Writer
{
	private:
		hid_t file_id;

		int open_file(std::string const &filename);
		hid_t dataset_create(hid_t &group_id, hid_t &dataspace_id, hsize_t count[2], std::string const &dataset);
		hid_t group_step_access(hid_t &file_id, long step);
		hid_t group_access(hid_t &file_id, std::string const &group);

	public:
		Writer(const std::string &filename);
		~Writer();

		int dump_serial(const Field_Data &field, long step);
		int dump_serial(std::string const &filename, long step, std::string const &group, const Field_Data &field);

		int write_attributes_parallel(const SimParams &simparams) const;

		template <class T>
		int writeattribute_file(std::string const &attr_name, T attr_value) const
		{
			writeattribute(file_id, attr_name, attr_value);
		}

		template <class T>
		int writeattribute(hid_t loc_id, std::string const &attr_name, T attr_value) const
		{
			hid_t type_id = H5T_NATIVE_DOUBLE;
			const std::type_info &type = typeid(attr_value);
			if (typeid(attr_value) == typeid(double)) {
				type_id = H5T_NATIVE_DOUBLE;
			} else if (typeid(attr_value) == typeid(int)) {
				type_id = H5T_NATIVE_INT;
			} else if (typeid(attr_value) == typeid(long)) {
				type_id = H5T_NATIVE_LONG;
			} else {
				printf("Couldn't write attribute\n");
			}

			hid_t attr_id, dataspace_id;
			herr_t status;
			const hsize_t mysize = 1;

			// ==================================
			// Create an appropriate dataspace
			// ==================================
			dataspace_id = H5Screate(H5S_SIMPLE);
			status = H5Sset_extent_simple(dataspace_id, 1, &mysize, &mysize);

			// ==================================
			// Create and write attribute
			// ==================================
			attr_id = H5Acreate(loc_id, attr_name.c_str(), type_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			H5Sclose(dataspace_id);
			status = H5Awrite(attr_id, type_id, &attr_value);
			H5Aclose(attr_id);

			return 0;
		}
};

#endif
