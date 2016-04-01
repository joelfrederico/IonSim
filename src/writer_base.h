#ifndef __WRITER_BASE_H_INCLUDED__
#define __WRITER_BASE_H_INCLUDED__

#include <string>
#include <hdf5.h>
#include "field_data.h"
#include "simparams.h"

class WriterBase
{
	protected:
		std::string _filename;
		hid_t file_id;

		hid_t dataset_create(hid_t &group_id, hid_t &dataspace_id, int rank, hsize_t *count, std::string const &dataset);
		hid_t group_step_access(hid_t &file_id, long step);
		hid_t group_access(hid_t &loc_id, std::string const &group);
		hid_t group_access(hid_t &loc_id, const char *group);

	public:
		WriterBase(const std::string &filename);
		~WriterBase();

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
				printf("Not a valid attribute type.\n");
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
			/* attr_id = -1; */
			if (attr_id < 0) {
				std::cout << "Attribute not created!" << std::endl;
				return -1;
			} else {
				H5Sclose(dataspace_id);
				status = H5Awrite(attr_id, type_id, &attr_value);
				H5Aclose(attr_id);
			}

			return 0;
		}
};

#endif
