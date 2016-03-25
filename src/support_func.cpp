#include "consts.h"
#include "fields.h"
#include "support_func.h"
#include <hdf5.h>
#include <iomanip>
#include <math.h>
#include <mpi.h>
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

	hid_t open_file(std::string const &filename)
	{
		hid_t file_id;
		return H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	}

	hid_t group_step_access(hid_t &file_id, long step)
	{
		char buf[10];
		sprintf(buf, "Step_%04ld", step);
		return group_access(file_id, buf);
	}

	hid_t group_access(hid_t &file_id, std::string const &group)
	{
		return group_access(file_id, group.c_str());
	}

	hid_t group_access(hid_t &file_id, const char *group)
	{
		herr_t status;
		hid_t group_id;
		H5G_info_t objinfo;

		// ==================================
		// Access or create a new group
		// ==================================
		H5Eset_auto(NULL, NULL, NULL);
		group_id = H5Gopen(file_id, group, H5P_DEFAULT);
		if (group_id < 0)
		{
			group_id = H5Gcreate(file_id, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Gget_info_by_name(file_id, group, &objinfo, H5P_DEFAULT);
		}
		return group_id;
	}

	hid_t dataset_create(hid_t &group_id, hid_t &dataspace_id, hsize_t count[2], std::string const &dataset)
	{
		hid_t dataset_id;

		dataspace_id = H5Screate_simple(2, count, NULL);
	
		// ==================================
		// Create dataset
		// ==================================
		dataset_id = H5Dcreate(group_id, dataset.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
		return dataset_id;
	}

	int dump_serial(std::string const &filename, long step, std::string const &group, const Field &field)
	{
		// ==================================
		// Initialize all variables
		// ==================================
		double *xbuf;
		double *ybuf;
		std::string x_dataset, y_dataset;
		long x_len = field.x_pts;
		long y_len = field.y_pts;

		int n_write = MAX_N_WRITE;

		hid_t file_id;
		hid_t step_group_id, group_id;
		hid_t x_dataspace_id, y_dataspace_id;
		hid_t x_dataset_id, y_dataset_id;
		hid_t x_memspace_id, y_memspace_id;
		hsize_t count[2];
		hsize_t offset[2];
		std::string temp_str;
		herr_t status;
	
		// ==================================
		// Open/create file
		// ==================================
		file_id = open_file(filename);
	
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
		count[0] = x_len;
		count[1] = y_len;
		x_dataset = "Ex";
		x_dataset_id = dataset_create(group_id, x_dataspace_id, count, x_dataset); // Updates dataspace_id
		x_memspace_id = H5Screate_simple(2, count, NULL);

		count[0] = x_len;
		count[1] = y_len;
		y_dataset = "Ey";
		y_dataset_id = dataset_create(group_id, y_dataspace_id, count, y_dataset); // Updates dataspace_id
		y_memspace_id = H5Screate_simple(2, count, NULL);

		// ==================================
		// Write dataset rows
		// ==================================
		xbuf = new double[x_len*y_len];
		ybuf = new double[x_len*y_len];
		for (int i=0; i < x_len; i++)
		{
			for (int j=0; j < y_len; j++)
			{
				xbuf[i + j*x_len] = field.Ex_ind(i, j);
				ybuf[i + j*x_len] = field.Ey_ind(i, j);
			}
		}
		H5Dwrite(x_dataset_id, H5T_NATIVE_DOUBLE, x_memspace_id, x_dataspace_id, H5P_DEFAULT, xbuf);
		H5Dwrite(y_dataset_id, H5T_NATIVE_DOUBLE, y_memspace_id, y_dataspace_id, H5P_DEFAULT, ybuf);

		delete [] xbuf;
		delete [] ybuf;
	
		// ==================================
		// Close file
		// ==================================
		/* hdf5_cleanup(x_dataspace_id, y_dataspace_id, x_dataset_id, y_dataset_id, file_id); */
		H5Gclose(step_group_id);
		H5Gclose(group_id);
		H5Sclose(x_dataspace_id);
		H5Sclose(y_dataspace_id);
		H5Dclose(x_dataset_id);
		H5Dclose(y_dataset_id);
		H5Sclose(x_memspace_id);
		H5Sclose(y_memspace_id);
		H5Fclose(file_id);
	
		return 0;
	}

	int dump_parallel(std::string const &filename, long step, std::string const &dataset, MPI::Intracomm &comm, const Parts &parts)
	{
		// ==================================
		// Initialize all variables
		// ==================================
		MPI::Info info;
		int p  = comm.Get_size();
		int id = comm.Get_rank();

		long n_pts             = parts.n_pts;
		const double_vec * _x  = &parts.x;
		const double_vec * _xp = &parts.xp;
		const double_vec * _y  = &parts.y;
		const double_vec * _yp = &parts.yp;
		const double_vec * _z  = &parts.z;
		const double_vec * _zp = &parts.zp;
		double *buf;

		int n_write = MAX_N_WRITE;

		hid_t plist_dx_id;
		hid_t file_id;
		hid_t step_group_id;
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
		step_group_id    = group_step_access(file_id, step);

		// ==================================
		// Create dataset collectively
		// ==================================
		// The total array of particles is (n_pts*num_processes, 6) in size
		// Remember that each slave has n_pts of different particles!
		count[0] = n_pts*p;
		count[1] = 6;
		dataset_id = dataset_create(step_group_id, dataspace_id, count, dataset);

		// ==================================
		// Write dataset rows
		// ==================================
		plist_dx_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_dx_id, H5FD_MPIO_COLLECTIVE);
	
		int i=0;
		if (n_pts < n_write) n_write = n_pts;
		/* double buf[n_write*6]; */
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
		delete [] buf;
	
		// ==================================
		// Close file
		// ==================================
		H5Pclose(plist_dx_id);
		H5Gclose(step_group_id);
		H5Sclose(dataspace_id);
		H5Dclose(dataset_id);
		H5Sclose(memspace_id);
		H5Fclose(file_id);

		return 0;
	}

	int dump_serial(std::string const &filename, long step, std::string const &group, std::string const &dataset, const Parts &ebeam)
	{
		// ==================================
		// Initialize all variables
		// ==================================
		long n_pts             = ebeam.n_pts;
		const double_vec * _x  = &ebeam.x;
		const double_vec * _xp = &ebeam.xp;
		const double_vec * _y  = &ebeam.y;
		const double_vec * _yp = &ebeam.yp;
		const double_vec * _z  = &ebeam.z;
		const double_vec * _zp = &ebeam.zp;
		double *buf;

		int n_write = MAX_N_WRITE;

		hid_t step_group_id, group_id, file_id, dataspace_id;
		hid_t dataset_id;
		hid_t memspace_id;
		hsize_t count[2];
	
		herr_t status;
	
		// ==================================
		// Open/create file
		// ==================================
		file_id = open_file(filename);
	
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
		count[0] = n_pts;
		count[1] = 6;
		dataset_id = dataset_create(group_id, dataspace_id, count, dataset);

	
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

		H5Fclose(file_id);
	
		return 0;
	}

	int overwrite_file_parallel(std::string const &filename, MPI::Intracomm &slave_comm_id)
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

	int overwrite_file_serial(std::string const &filename)
	{
		hid_t file_id;
		MPI::Info info;

		// ==================================
		// Create a new file collectively
		// ==================================
		file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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

	int sendloop(const int &message)
	{
		int buf;
		buf = message;
		MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
		return 0;
	}

	int loop_get_fields(Field &field)
	{
		sendloop(LOOP_GET_EFIELD);
		field.recv_field_others();
		return 0;
	}

	int loop_push_ions(Field &field)
	{
		sendloop(LOOP_PUSH_IONS);
		int p = MPI::COMM_WORLD.Get_size();
		for (int id=1; id < p; id++) {
			field.send_field(id);
		}
		return 0;
	}

	int sendloop(const int &message, int step)
	{
		sendloop(message);
		MPI::COMM_WORLD.Bcast(&step, 1, MPI::INT, 0);
		return 0;
	}
}
