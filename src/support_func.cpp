#include "support_func.h"
#include "consts.h"
#include <hdf5.h>
#include <mpi.h>
#include "fields.h"
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>

// Number of entries to write at once
// because HDF5 crashes if you try to
// do too much all at once.
const int MAX_N_WRITE = 1e5;

namespace ionsim
{
	// ==================================
	// Methods
	// ==================================
	double GeV2gamma(double GeV)
	{
		return GeV * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT / ELECTRON_REST_ENERGY;
	}

	double gamma2GeV(double gamma)
	{
		return gamma * ELECTRON_REST_ENERGY / (1e9 * GSL_CONST_MKSA_ELECTRON_VOLT);
	}

	double gaussian()
	{
		return 0;
	}

	hid_t open_file_parallel(std::string const &filename, MPI::Intracomm &slave_comm_id)
	{
		// ==================================
		// Set up file access property list
		// ==================================
		MPI::Info info;
		hid_t file_id;
		hid_t plist_file_id;

		plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_file_id, slave_comm_id, info);
	
		// ==================================
		// Open the file
		// ==================================
		file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_file_id);
		H5Pclose(plist_file_id);
		return file_id;
	}

	hid_t open_file(std::string const &filename, MPI::Intracomm &slave_comm_id)
	{
		hid_t file_id;
		return H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	}

	hid_t group_access(hid_t &file_id, std::string const &group)
	{
		herr_t status;
		hid_t group_id;

		// ==================================
		// Access or create a new group
		// ==================================
		H5Eset_auto(NULL, NULL, NULL);
		status = H5Gget_info_by_name(file_id, group.c_str(), NULL, H5P_DEFAULT);
		if (status >= 0)
		{
			group_id = H5Gopen(file_id, group.c_str(), H5P_DEFAULT);
		} else {
			group_id = H5Gcreate(file_id, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		return group_id;
	}

	hid_t dataset_create(hid_t &group_id, hid_t &dataspace_id, hsize_t count[2], std::string const &dataset)
	{
		std::stringstream dataset_str;
		hid_t dataset_id;

		dataspace_id = H5Screate_simple(2, count, NULL);
	
		// ==================================
		// Create dataset
		// ==================================
		dataset_str << "step_" << std::setfill('0') << std::setw(3) << dataset;
		dataset_id = H5Dcreate(group_id, dataset_str.str().c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Gclose(group_id);
	
		return dataset_id;
	}

	int hdf5_cleanup(hid_t &dataspace_id, hid_t &dataset_id, hid_t &file_id)
	{
		H5Sclose(dataspace_id);
		H5Dclose(dataset_id);
		H5Fclose(file_id);
		return 0;
	}

	int serial_dump(std::string const &filename, long step, std::string const &group, std::string const &dataset, MPI::Intracomm &comm, const Field &field)
	{
		// ==================================
		// Initialize all variables
		// ==================================
		MPI::Info info;
		int p  = comm.Get_size();
		int id = comm.Get_rank();

		/* long n_pts             = ebeam.n_pts; */
		/* const double_vec * _x  = &ebeam.x; */
		/* const double_vec * _xp = &ebeam.xp; */
		/* const double_vec * _y  = &ebeam.y; */
		/* const double_vec * _yp = &ebeam.yp; */
		/* const double_vec * _z  = &ebeam.z; */
		/* const double_vec * _zp = &ebeam.zp; */

		long x_len = field.x_pts;
		long y_len = field.y_pts;

		int n_write = MAX_N_WRITE;

		hid_t plist_dx_id;
		hid_t file_id;
		hid_t step_group_id, group_id, group_field_x_id, group_field_y_id;
		hid_t dataspace_id;
		hid_t dataset_id;
		hid_t memspace_id;
		hsize_t count[2] = {x_len, y_len};
		hsize_t offset[2];
		std::string temp_str;
		herr_t status;
	
		// ==================================
		// Open/create file
		// ==================================
		file_id = open_file(filename, comm);
	
		// ==================================
		// Access or create a new group
		// ==================================
		temp_str         = "step_" + std::to_string(step);
		step_group_id    = group_access(file_id, temp_str);

		group_id         = group_access(file_id, group);

		temp_str         = "field_x";
		group_field_x_id = group_access(group_id, temp_str);

		// ==================================
		// Create dataset collectively
		// ==================================
		// The total array of particles is (n_pts*num_processes, 6) in size
		// Remember that each slave has n_pts of different particles!
		count[0] = x_len;
		count[1] = y_len;
		dataset_id = dataset_create(group_field_x_id, dataspace_id, count, dataset); // Updates dataspace_id

		// ==================================
		// Write dataset rows
		// ==================================
		memspace_id = H5Screate_simple(2, count, NULL);
		double buf[x_len*y_len];
		for (int i=0; i < x_len; i++)
		{
			for (int j=0; j < y_len; j++)
			{
				buf[i + j*x_len] = field.Ex_ind(i, j);
			}
		}
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, buf);
	
		// ==================================
		// Close file
		// ==================================
		hdf5_cleanup(dataspace_id, dataset_id, file_id);
	
		return 0;
	}

	int dump(std::string const &filename, std::string const &group, std::string const &dataset, MPI::Intracomm &comm, const Parts &ebeam)
	{
		// ==================================
		// Initialize all variables
		// ==================================
		MPI::Info info;
		int p  = comm.Get_size();
		int id = comm.Get_rank();

		long n_pts             = ebeam.n_pts;
		const double_vec * _x  = &ebeam.x;
		const double_vec * _xp = &ebeam.xp;
		const double_vec * _y  = &ebeam.y;
		const double_vec * _yp = &ebeam.yp;
		const double_vec * _z  = &ebeam.z;
		const double_vec * _zp = &ebeam.zp;

		int n_write = MAX_N_WRITE;

		hid_t plist_dx_id;
		hid_t file_id;
		hid_t group_id;
		hid_t dataspace_id;
		hid_t dataset_id;
		hid_t memspace_id;
		hsize_t count[2] = {n_pts*p, 6};
		hsize_t offset[2];
	
		herr_t status;
	
		// ==================================
		// Open/create file
		// ==================================
		file_id = open_file_parallel(filename, comm);
	
		// ==================================
		// Access or create a new group
		// ==================================
		group_id = group_access(file_id, group);

		// ==================================
		// Create dataset collectively
		// ==================================
		// The total array of particles is (n_pts*num_processes, 6) in size
		// Remember that each slave has n_pts of different particles!
		count[0] = n_pts*p;
		count[1] = 6;
		dataset_id = dataset_create(group_id, dataspace_id, count, dataset);

		/* hsize_t offset_init[2] = {n_pts, 0}; */
	
		// ==================================
		// Write dataset rows
		// ==================================
		plist_dx_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_dx_id, H5FD_MPIO_COLLECTIVE);
	
		int i=0;
		if (n_pts < n_write) n_write = n_pts;
		double buf[n_write*6];
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
			offset[0] = (id * n_pts + i);
			offset[1] = 0;

			memspace_id = H5Screate_simple(2, count, NULL);
	
			H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	
	
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
	
			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, plist_dx_id, buf);
		}
	
		// ==================================
		// Close file
		// ==================================
		H5Pclose(plist_dx_id);
		hdf5_cleanup(dataspace_id, dataset_id, file_id);
	
		return 0;
	}

	int overwrite_file(std::string const &filename, MPI::Intracomm &slave_comm_id)
	{
		hid_t file_id;
		hid_t plist_file_id;
		MPI::Info info;

		// ==================================
		// Set up file access property list
		// ==================================
		plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_file_id, slave_comm_id, info);

		// ==================================
		// Create a new file collectively
		// ==================================
		file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file_id);
		H5Pclose(plist_file_id);
		H5Fclose(file_id);
		
		return 0;
	}

	double ** alloc_2d_array(long rowCount, long colCount)
	{
		double ** out = new double*[rowCount];
		for (long i=0; i < rowCount; i++)
		{
			out[i] = new double[colCount];
		}
		return out;
	}

	int dealloc_2d_array(double ** (&arr), long rowCount)
	{
		for (long i=0; i < rowCount; i++)
		{
			delete [] arr[i];
		}
		delete [] arr;
		return 0;
	}

	int sendloop(const int * message)
	{
		int buf;
		buf = *message;
		MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
		return 0;
	}

	int sendloop(const int * message, int step)
	{
		int buf;
		sendloop(message);
		MPI::COMM_WORLD.Bcast(&step, 1, MPI::INT, 0);
		return 0;
	}
}
