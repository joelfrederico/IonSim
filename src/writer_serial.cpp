#include "writer_serial.h"
#include "parts.h"
#include "support_func.h"
#include "field_data.h"
#include <ionsim.h>

WriterSerial::WriterSerial(const std::string &filename) :
	WriterBase(filename)
{
	_init(filename, false);
}

WriterSerial::WriterSerial(const std::string &filename, bool overwrite) :
	WriterBase(filename)
{
	_init(filename, overwrite);
}

int WriterSerial::_init(const std::string &filename, bool overwrite)
{
	if (overwrite)
	{
		overwrite_file_serial();
	} else {
		open_file();
	}

	return 0;
}

int WriterSerial::open_file()
{
	file_id = H5Fopen(_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	return 0;
}

int WriterSerial::overwrite_file_serial()
{
	// ==================================
	// Create a new file
	// ==================================
	file_id = H5Fcreate(_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	return 0;
}

int WriterSerial::writedata(long step, const std::string &group_str, const std::string &dataset_str, const Parts &parts)
{
	// ==================================
	// Initialize all variables
	// ==================================
	long long n_pts             = parts.n_pts;
	const ldouble_vec * _x  = &parts.x;
	const ldouble_vec * _xp = &parts.xp;
	const ldouble_vec * _y  = &parts.y;
	const ldouble_vec * _yp = &parts.yp;
	const ldouble_vec * _z  = &parts.z;
	const ldouble_vec * _zp = &parts.zp;
	double *buf;

	int n_write = MAX_N_WRITE;

	DataspaceCreate *memspace;
	hsize_t count[2];
	herr_t status;

	// ==================================
	// Access or create a new group
	// ==================================
	GroupStepAccess step_group = GroupStepAccess(file_id, step);

	GroupAccess group = GroupAccess(step_group.group_id, group_str);

	// ==================================
	// Create dataset collectively
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	int rank = 2;
	count[0] = n_pts;
	count[1] = 6;
	DatasetAccess dataset(group.group_id, dataset_str, rank, count);

	// ==================================
	// Write dataset rows
	// ==================================
	int i=0;
	if (n_pts < n_write) n_write = n_pts;
	buf = new double[n_write*6];
	for (int j=0; i < n_pts; j++)
	{
		if (n_pts - i < n_write)
		{
			n_write = n_pts - i;
		}
		// ==================================
		// Write to hyperslab
		// ==================================
		count[0]  = n_write;
		count[1]  = 6;

		memspace = new DataspaceCreate(2, count);

		for (int k=0; k < n_write; k++)
		{
			buf[k*6]     = (*_x) [i];
			buf[k*6 + 1] = (*_xp)[i];
			buf[k*6 + 2] = (*_y) [i];
			buf[k*6 + 3] = (*_yp)[i];
			buf[k*6 + 4] = (*_z) [i];
			buf[k*6 + 5] = (*_zp)[i];
			i++;
		}

		status = H5Dwrite(dataset.dataset_id, H5T_NATIVE_DOUBLE, (*memspace).dataspace_id, dataset.dataspace_id, H5P_DEFAULT, buf);

		delete memspace;
	}

	// ==================================
	// Close file
	// ==================================
	delete [] buf;

	return 0;
}

hid_t _write_grid(hid_t loc_id, std::string attr_name, double *attr_array, const hsize_t size)
{
	herr_t status;
	hid_t attr_id;

	DataspaceCreate x_grid_space(H5S_SCALAR);
	status = H5Sset_extent_simple(x_grid_space.dataspace_id, 1, &size, NULL);

	attr_id = H5Acreate(loc_id, attr_name.c_str(), H5T_NATIVE_DOUBLE, x_grid_space.dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, attr_array);

	H5Aclose(attr_id);

	return status;
}

int write_data(GroupAccess *group, Field_Data *field, double *buf, std::string dataset_str)
{
	// ==================================
	// Initialize all variables
	// ==================================
	double *_buf;
	hsize_t count[3];
	hsize_t offset[3];
	std::string temp_str;
	herr_t status;

	// ==================================
	// Create dataset
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	int rank = 3;
	int x_len = field->x_pts;
	int y_len = field->y_pts;
	int z_len = field->z_pts;

	/* std::cout << "x_len: " << x_len << std::endl; */
	/* std::cout << "y_len: " << y_len << std::endl; */
	/* std::cout << "z_len: " << z_len << std::endl; */

	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;
	DatasetAccess dataset(group->group_id, dataset_str, rank, count); // Updates dataspace_id
	DataspaceCreate memspace(rank, count);

	// ==================================
	// Write dataset rows
	// ==================================
	_buf = new double[x_len*y_len*z_len];

	for (int i=0; i < x_len; i++) {
		for (int j=0; j < y_len; j++) {
			for (int k=0; k < z_len; k++) {
				_buf[k + z_len*(i + j*x_len)] = buf[field->_index(i, j, k)];
			}
		}
	}
	status = H5Dwrite(dataset.dataset_id, H5T_NATIVE_DOUBLE, memspace.dataspace_id, dataset.dataspace_id, H5P_DEFAULT, _buf);

	if (status > 0) std::cout << "Something bad" << std::endl;

	delete [] _buf;

	return 0;
}

int WriterSerial::writedata(long step, Field_Data &field)
{
	// ==================================
	// Initialize all variables
	// ==================================
	int x_len = field.x_pts;
	int y_len = field.y_pts;
	int z_len = field.z_pts;
	const std::string group_str = "field"; 

	// ==================================
	// Access or create a new group
	// ==================================
	GroupStepAccess step_group = GroupStepAccess(file_id, step);
	GroupAccess group = GroupAccess(step_group.group_id, group_str);

	// ==================================
	// Write field attributes
	// ==================================
	_write_grid(group.group_id, "x_grid", field.x_grid, field.x_pts);
	_write_grid(group.group_id, "y_grid", field.y_grid, field.y_pts);
	_write_grid(group.group_id, "z_grid", field.z_grid, field.z_pts);

	AttributeCreate z_pts(group.group_id, "z_pts"          , field.z_pts);

	// ==================================
	// Write field
	// ==================================
	write_data(&group, &field, field.Ex_data, "Ex");
	write_data(&group, &field, field.Ey_data, "Ey");
	write_data(&group, &field, field.Ez_data, "Ez");

	write_data(&group, &field, field.Bx_data, "Bx");
	write_data(&group, &field, field.By_data, "By");
	write_data(&group, &field, field.Bz_data, "Bz");

	return 0;
}
