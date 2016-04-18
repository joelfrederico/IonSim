#include <config.h>
#include <stdio.h>
#include <hdf5.h>
#include <iomanip>
#include <argp.h>
#include <string>

const char *argp_program_version = "v0.1";
static char args_doc[] = "FILE";

herr_t print_attrs(hid_t loc_id, const char *name, const H5A_info_t *info, void *opdata)
{
	hid_t attr_id;
	hid_t attr_type;
	herr_t status;
	int i_buf;
	double d_buf;

	attr_id   = H5Aopen(loc_id, name, H5P_DEFAULT);
	attr_type = H5Aget_type(attr_id);

	switch(H5Tget_class(attr_type))
	{
		case H5T_INTEGER:
			status = H5Aread(attr_id, attr_type, &i_buf);
			std::cout << std::left << std::setw(15) << name << ": " << i_buf << std::endl;
			break;
		case H5T_FLOAT:
			status = H5Aread(attr_id, attr_type, &d_buf);
			std::cout << std::left << std::setw(15) << name << ": " << d_buf << std::endl;
			break;
		default:
			std::cout << std::left << std::setw(15) << name << ": " << "(attribute type not known)" << std::endl;
	}

	H5Tclose(attr_type);
	H5Aclose(attr_id);
	return 0;
}

struct arguments
{
	std::string file;
};

error_t myparser(int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments;
	arguments = (struct arguments*)state->input;

	switch (key)
	{
		case ARGP_KEY_ARG:
			arguments->file = arg;
			break;
		default:
			return ARGP_ERR_UNKNOWN;
	}

	return 0;
}

int main(int argc, char **argv)
{
	int arg_err;

	struct arguments args;
	args.file = "output.h5";

	struct argp argp = {0, myparser, args_doc, 0};
	argp_parse(&argp, argc, argv, 0, 0, &args);

	if (arg_err != 0)
	{
		std::cout << "Error: " << arg_err << std::endl;
		return arg_err;
	}

	std::cout << "Checking file: " << args.file << std::endl;
	std::cout << "=========================================" << std::endl;

	hid_t file_id, root_obj;
	herr_t status;
	H5O_info_t root_info;
	hsize_t n = 0;

	file_id = H5Fopen(args.file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	status = H5Oget_info(file_id, &root_info);

	std::cout << "Number of Attrs: " << root_info.num_attrs << std::endl;
	std::cout << "-----------------------------------------" << std::endl;
	H5Aiterate(file_id, H5_INDEX_NAME, H5_ITER_INC, &n, print_attrs, NULL);

	H5Fclose(file_id);

	return 0;
}
