#include "support_func.h"
#include "consts.h"
#include "hdf5.h"
#include "mpi.h"
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>

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

	hid_t open_file(std::string const &filename, MPI::Intracomm &slave_comm_id)
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

	int dump(std::string const &filename, std::string const &group, std::string const &dataset, MPI::Intracomm &comm, const Parts &ebeam)
	{
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

		int n_write = 1e5;

		hid_t plist_dx_id;
		hid_t file_id;
		hid_t group_id;
		hid_t dataspace_id;
		hid_t dataset_id;
		std::stringstream dataset_str;
		hid_t memspace_id;
		hsize_t count[2] = {n_pts*p, 6};
		hsize_t offset[2];
	
		herr_t status;
	
		file_id = open_file(filename, comm);
	
		// ==================================
		// Access or create a new group
		// ==================================
		H5Eset_auto(NULL, NULL, NULL);
		status = H5Gget_objinfo(file_id, group.c_str(), 0, NULL);
		std::stringstream group_path;
		group_path << "/" << group;
		if (status == 0)
		{
			group_id = H5Gopen(file_id, group.c_str(), H5P_DEFAULT);
		} else {
			group_id = H5Gcreate(file_id, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
	
		// ==================================
		// Create dataspace collectively
		// ==================================
		/* memcpy(count, (hsize_t [2]){n_pts*comm.Get_size(), 6}, 2*sizeof(hsize_t)); */
		count[0] = n_pts*p;
		count[1] = 6;
		dataspace_id = H5Screate_simple(2, count, NULL);
	
		// ==================================
		// Create dataset
		// ==================================
		dataset_str << "step_" << std::setfill('0') << std::setw(3) << dataset;
		dataset_id = H5Dcreate(group_id, dataset_str.str().c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Gclose(group_id);
	
		hsize_t offset_init[2] = {n_pts, 0};
	
		plist_dx_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_dx_id, H5FD_MPIO_COLLECTIVE);
	
		int i=0;
		if (n_pts < n_write) n_write = n_pts;
		double buf[n_write*6];
		for (int j=0; i < n_pts; j++)
		{
			/* if (id == 0) printf("On chunk %d...\n", j+1); */
	
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
		H5Sclose(dataspace_id);
		H5Dclose(dataset_id);
		H5Fclose(file_id);
	
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


}
