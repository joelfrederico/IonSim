#include <gsl/gsl_const_mksa.h>
#include <math.h>
#include <string>
#include "hdf5.h"
#include "classes.h"

namespace ionsim
{
	const double ELECTRON_REST_ENERGY = GSL_CONST_MKSA_MASS_ELECTRON * pow(GSL_CONST_MKSA_SPEED_OF_LIGHT, 2);

	double gamma2GeV(double gamma)
	{
		return gamma * ELECTRON_REST_ENERGY / (1e9 * GSL_CONST_MKSA_ELECTRON_VOLT);
	}

	double GeV2gamma(double GeV)
	{
		return GeV * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT / ELECTRON_REST_ENERGY;
	}

	int dump(std::string const &filename, MPI::Intracomm &comm, Parts *ebeam)
	{
		long n_pts             = (*ebeam).n_pts();
		const double_vec * _x  = (*ebeam).x();
		const double_vec * _xp = (*ebeam).xp();
		const double_vec * _y  = (*ebeam).y();
		const double_vec * _yp = (*ebeam).yp();
		const double_vec * _z  = (*ebeam).z();
		const double_vec * _zp = (*ebeam).zp();

		hid_t plist_file_id;
		hid_t plist_dx_id;
		hid_t file_id;
		hid_t group_id;
		hid_t dataspace_id;
		hid_t dataset_id;
		hid_t memspace_id;
		hsize_t count[2];
		hsize_t offset[2];
	
		herr_t status;
	
		MPI::Info info;
		int p  = comm.Get_size();
		int id = comm.Get_rank();
	
		// ==================================
		// Set up file access property list
		// ==================================
		plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_file_id, comm, info);
	
		// ==================================
		// Create a new file collectively,
		// ==================================
		file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file_id);
		H5Pclose(plist_file_id);
	
		// ==================================
		// Create group collectively
		// ==================================
		group_id = H5Gcreate(file_id, "ebeam", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
		// ==================================
		// Create dataspace collectively
		// ==================================
		memcpy(count, (hsize_t [2]){n_pts*comm.Get_size(), 6}, 2*sizeof(hsize_t));
		dataspace_id = H5Screate_simple(2, count, NULL);
	
		// ==================================
		// Create dataset
		// ==================================
		dataset_id = H5Dcreate(group_id, "particles", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Gclose(group_id);
	
		hsize_t offset_init[2] = {n_pts, 0};
	
		plist_dx_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_dx_id, H5FD_MPIO_COLLECTIVE);
	
		int n_write=1e5;
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
			memcpy(count, (hsize_t [2]){n_write, 6}, 2*sizeof(hsize_t));
			memcpy(offset, (hsize_t [2]){id * n_pts + i, 0}, 2*sizeof(hsize_t));
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
}
