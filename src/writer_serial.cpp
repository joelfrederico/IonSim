#include "writer_serial.h"
#include "parts.h"
#include "support_func.h"
#include "field_data.h"
#include "hdf5_classes.h"

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
	long n_pts             = parts.n_pts;
	const double_vec * _x  = &parts.x;
	const double_vec * _xp = &parts.xp;
	const double_vec * _y  = &parts.y;
	const double_vec * _yp = &parts.yp;
	const double_vec * _z  = &parts.z;
	const double_vec * _zp = &parts.zp;
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

int WriterSerial::writedata(long step, Field_Data &field)
{
	// ==================================
	// Initialize all variables
	// ==================================
	double *xbuf;
	double *ybuf;
	double *zbuf;
	std::string x_dataset_str, y_dataset_str, z_dataset_str;
	long x_len = field.x_pts;
	long y_len = field.y_pts;
	long z_len = field.z_pts;
	const std::string group_str = "field"; 

	int n_write = MAX_N_WRITE;

	hsize_t count[3];
	hsize_t offset[3];
	std::string temp_str;
	herr_t status;


	// ==================================
	// Access or create a new group
	// ==================================
	GroupStepAccess step_group = GroupStepAccess(file_id, step);

	GroupAccess group = GroupAccess(step_group.group_id, group_str);

	// ==================================
	// Create dataset
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	int rank = 3;
	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;
	x_dataset_str = "Ex";
	DatasetAccess x_dataset(group.group_id, x_dataset_str, rank, count); // Updates dataspace_id
	DataspaceCreate x_memspace(rank, count);

	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;
	y_dataset_str = "Ey";
	DatasetAccess y_dataset(group.group_id, y_dataset_str, rank, count); // Updates dataspace_id
	DataspaceCreate y_memspace(rank, count);

	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;
	z_dataset_str = "Ez";
	DatasetAccess z_dataset(group.group_id, z_dataset_str, rank, count); // Updates dataspace_id
	DataspaceCreate z_memspace(rank, count);

	// ==================================
	// Write dataset rows
	// ==================================
	xbuf = new double[x_len*y_len*z_len];
	ybuf = new double[x_len*y_len*z_len];
	zbuf = new double[x_len*y_len*z_len];

	for (int i=0; i < x_len; i++) {
		for (int j=0; j < y_len; j++) {
			for (int k=0; k < z_len; k++) {
				xbuf[k + z_len*(i + j*x_len)] = field.Ex_ind(i, j, k);
				ybuf[k + z_len*(i + j*x_len)] = field.Ey_ind(i, j, k);
				zbuf[k + z_len*(i + j*x_len)] = field.Ez_ind(i, j, k);
			}
		}
	}
	status = H5Dwrite(x_dataset.dataset_id, H5T_NATIVE_DOUBLE, x_memspace.dataspace_id, x_dataset.dataspace_id, H5P_DEFAULT, xbuf);
	status = H5Dwrite(y_dataset.dataset_id, H5T_NATIVE_DOUBLE, y_memspace.dataspace_id, y_dataset.dataspace_id, H5P_DEFAULT, ybuf);
	status = H5Dwrite(z_dataset.dataset_id, H5T_NATIVE_DOUBLE, z_memspace.dataspace_id, z_dataset.dataspace_id, H5P_DEFAULT, zbuf);

	if (status > 0) std::cout << "Something bad" << std::endl;

	delete [] xbuf;
	delete [] ybuf;
	delete [] zbuf;

	return 0;
}
