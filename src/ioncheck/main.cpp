#include <config.h>
#include <stdio.h>
#include <hdf5.h>
#include <iomanip>
/* #include <argp.h> */
#include <string>
#include <stdlib.h>

const char *argp_program_version = "v0.1";
static char args_doc[] = "FILE";
const int WIDTH = 17;

herr_t print_attrs(hid_t loc_id, const char *name, const H5A_info_t *info, void *opdata)
{
	hid_t attr_id, attr_type;
	herr_t status;
	int i_buf;
	double d_buf;
	char *c_buf;

	attr_id   = H5Aopen(loc_id, name, H5P_DEFAULT);
	attr_type = H5Aget_type(attr_id);

	switch(H5Tget_class(attr_type))
	{
		case H5T_INTEGER:
			status = H5Aread(attr_id, attr_type, &i_buf);
			std::cout << std::left << std::setw(WIDTH) << name << ": " << i_buf << std::endl;
			break;

		case H5T_FLOAT:
			status = H5Aread(attr_id, attr_type, &d_buf);
			std::cout << std::left << std::setw(WIDTH) << name << ": " << d_buf << std::endl;
			break;

		case H5T_STRING:
			status = H5Aread(attr_id, attr_type, &c_buf);

			std::cout << std::left << std::setw(WIDTH) << name << ": \"" << c_buf << "\"" << std::endl;

			/* delete c_buf; */
			free(c_buf);

			break;

		default:
			std::cout << std::left << std::setw(WIDTH) << name << ": " << "(attribute type not known)" << std::endl;
	}

	H5Tclose(attr_type);
	H5Aclose(attr_id);
	return 0;
}

int main(int argc, char **argv)
{
	const char *file;
	char *version;

	// ==================================
	// Determine whether to open default
	// ==================================
	if (argc == 2)
	{
		file = argv[1];
		std::cout << std::left << std::setw(25) << "Checking file: " << file << std::endl;
	} else {
		file = "output.h5";
		std::cout << std::left << std::setw(25) << "Checking default file: " << file << std::endl;
	}

	// ==================================
	// Try to open file
	// ==================================
	hid_t file_id, root_obj;
	herr_t status;
	H5O_info_t root_info;
	hsize_t n = 0;

	status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	file_id = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (file_id < 0)
	{
		std::cout << "File could not be opened!" << std::endl;
		return -1;
	}

	// ==================================
	// Print version
	// ==================================
	hid_t attr_id, attr_type;

	attr_id   = H5Aopen(file_id, "version", H5P_DEFAULT);
	attr_type = H5Aget_type(attr_id);
	status = H5Aread(attr_id, attr_type, &version);
	H5Aclose(attr_id);

	std::cout << std::left << std::setw(25) << "File version: " << version << std::endl;

	free(version);

	// ==================================
	// Write out attributes
	// ==================================
	std::cout << "=========================================" << std::endl;

	status = H5Oget_info(file_id, &root_info);

	std::cout << "Number of Attrs: " << root_info.num_attrs << std::endl;
	std::cout << "-----------------------------------------" << std::endl;
	H5Aiterate(file_id, H5_INDEX_NAME, H5_ITER_INC, &n, print_attrs, NULL);

	// ==================================
	// Close file
	// ==================================
	H5Fclose(file_id);

	return 0;
}
