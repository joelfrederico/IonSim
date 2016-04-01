#include "writer_serial.h"
#include "parts.h"
#include "support_func.h"
#include "field_data.h"

WriterSerial::WriterSerial(const std::string &filename) : WriterBase(filename)
{
	open_file(filename);
}

int WriterSerial::open_file(std::string const &filename)
{
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	return 0;
}

int WriterSerial::overwrite_file_serial(std::string const &filename)
{
	// ==================================
	// Create a new file
	// ==================================
	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	return 0;
}

int WriterSerial::writedata(long step, std::string const &group, std::string const &dataset, const Parts &parts)
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

	int n_write = ionsim::MAX_N_WRITE;

	hid_t step_group_id, group_id, file_id, dataspace_id;
	hid_t dataset_id;
	hid_t memspace_id;
	hsize_t count[2];

	herr_t status;

	// ==================================
	// Access or create a new group
	// ==================================
	step_group_id    = group_step_access(file_id, step);

	group_id         = group_access(step_group_id, group);

	// ==================================
	// Create dataset collectively
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	int rank = 2;
	count[0] = n_pts;
	count[1] = 6;
	dataset_id = dataset_create(group_id, dataspace_id, rank, count, dataset);

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

		memspace_id = H5Screate_simple(2, count, NULL);

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

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, buf);
	}

	// ==================================
	// Close file
	// ==================================
	delete [] buf;

	H5Sclose(step_group_id);
	H5Sclose(dataspace_id);
	H5Gclose(group_id);

	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);
	H5Sclose(memspace_id);

	return 0;
}

int WriterSerial::writedata(long step, const Field_Data &field)
	{
		std::cout << "Here" << std::endl;
		// ==================================
		// Initialize all variables
		// ==================================
		double *xbuf;
		double *ybuf;
		double *zbuf;
		std::string x_dataset, y_dataset, z_dataset;
		long x_len = field.x_pts;
		long y_len = field.y_pts;
		long z_len = field.z_pts;
		std::string const group = "field"; 

		int n_write = ionsim::MAX_N_WRITE;

		hid_t step_group_id, group_id;
		hid_t x_dataspace_id, y_dataspace_id, z_dataspace_id;
		hid_t x_dataset_id, y_dataset_id, z_dataset_id;
		hid_t x_memspace_id, y_memspace_id, z_memspace_id;
		hsize_t count[3];
		hsize_t offset[3];
		std::string temp_str;
		herr_t status;
	
		// ==================================
		// Access or create a new group
		// ==================================
		step_group_id    = group_step_access(file_id, step);

		group_id         = group_access(step_group_id, group);

		// ==================================
		// Create dataset
		// ==================================
		// The total array of particles is (n_pts*num_processes, 6) in size
		// Remember that each slave has n_pts of different particles!
		int rank = 3;
		count[0] = x_len;
		count[1] = y_len;
		count[2] = z_len;
		x_dataset = "Ex";
		x_dataset_id = dataset_create(group_id, x_dataspace_id, rank, count, x_dataset); // Updates dataspace_id
		x_memspace_id = H5Screate_simple(rank, count, NULL);

		count[0] = x_len;
		count[1] = y_len;
		count[2] = z_len;
		y_dataset = "Ey";
		y_dataset_id = dataset_create(group_id, y_dataspace_id, rank, count, y_dataset); // Updates dataspace_id
		y_memspace_id = H5Screate_simple(rank, count, NULL);

		count[0] = x_len;
		count[1] = y_len;
		count[2] = z_len;
		z_dataset = "Ez";
		z_dataset_id = dataset_create(group_id, z_dataspace_id, rank, count, z_dataset); // Updates dataspace_id
		z_memspace_id = H5Screate_simple(rank, count, NULL);

		// ==================================
		// Write dataset rows
		// ==================================
		xbuf = new double[x_len*y_len*z_len];
		ybuf = new double[x_len*y_len*z_len];
		zbuf = new double[x_len*y_len*z_len];

		for (int i=0; i < x_len; i++) {
			for (int j=0; j < y_len; j++) {
				for (int k=0; k < z_len; k++) {
					xbuf[i + x_len*(j + k*y_len)] = field.Ex_ind(i, j, k);
					ybuf[i + x_len*(j + k*y_len)] = field.Ey_ind(i, j, k);
					zbuf[i + x_len*(j + k*y_len)] = field.Ez_ind(i, j, k);
				}
			}
		}
		status = H5Dwrite(x_dataset_id, H5T_NATIVE_DOUBLE, x_memspace_id, x_dataspace_id, H5P_DEFAULT, xbuf);
		status = H5Dwrite(y_dataset_id, H5T_NATIVE_DOUBLE, y_memspace_id, y_dataspace_id, H5P_DEFAULT, ybuf);
		status = H5Dwrite(z_dataset_id, H5T_NATIVE_DOUBLE, z_memspace_id, z_dataspace_id, H5P_DEFAULT, zbuf);

		if (status > 0) std::cout << "Something bad" << std::endl;

		delete [] xbuf;
		delete [] ybuf;
		delete [] zbuf;
	
		// ==================================
		// Close file
		// ==================================
		H5Gclose(step_group_id);
		H5Gclose(group_id);
		H5Sclose(x_dataspace_id);
		H5Sclose(y_dataspace_id);
		H5Sclose(z_dataspace_id);
		H5Dclose(x_dataset_id);
		H5Dclose(y_dataset_id);
		H5Dclose(z_dataset_id);
		H5Sclose(x_memspace_id);
		H5Sclose(y_memspace_id);
		H5Sclose(z_memspace_id);
	
		return 0;
	}
