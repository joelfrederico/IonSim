#include "writer_parallel.h"
#include <mpi.h>
#include "support_func.h"

WriterParallel::WriterParallel(const std::string &filename, MPI::Intracomm *comm_id) : WriterBase(filename)
{
	_init(filename, comm_id, false);
}

WriterParallel::WriterParallel(const std::string &filename, MPI::Intracomm *comm_id, bool overwrite) : WriterBase(filename)
{
	_init(filename, comm_id, overwrite);
}

int WriterParallel::_init(const std::string &filename, MPI::Intracomm *comm_id, bool overwrite)
{
	comm_id_ptr = comm_id;
	p = MPI::COMM_WORLD.Get_size();
	id = MPI::COMM_WORLD.Get_rank();

	if (overwrite) {
		overwrite_file_parallel(filename);
	} else {
		open_file(filename);
	}

	return 0;
}

int WriterParallel::write_attributes(const SimParams &simparams) const
{
	writeattribute(file_id, "n_e"       , simparams.n_e      );
	writeattribute(file_id, "n_ions"    , simparams.n_ions   );
	writeattribute(file_id, "q_tot"     , simparams.q_tot    );
	writeattribute(file_id, "radius"    , simparams.radius   );
	writeattribute(file_id, "length"    , simparams.length   );
	writeattribute(file_id, "E"         , simparams.E        );
	writeattribute(file_id, "emit_n"    , simparams.emit_n   );
	writeattribute(file_id, "n_p_cgs"   , simparams.n_p_cgs  );
	writeattribute(file_id, "m_ion_amu" , simparams.m_ion_amu);
	writeattribute(file_id, "sz"        , simparams.sz       );
	writeattribute(file_id, "sdelta"    , simparams.sdelta   );
	writeattribute(file_id, "dt"        , simparams.dt       );
	writeattribute(file_id, "n_steps"   , simparams.n_steps  );
	writeattribute(file_id, "n_field_x" , simparams.n_field_x);
	writeattribute(file_id, "n_field_y" , simparams.n_field_y);
	writeattribute(file_id, "n_field_z" , simparams.n_field_z);

	return 0;
}

int WriterParallel::open_file(std::string const &filename)
{
	// ==================================
	// Set up file access property list
	// ==================================
	MPI::Info info;
	hid_t plist_file_id;

	plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_file_id, *comm_id_ptr, info);

	// ==================================
	// Open the file
	// ==================================
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_file_id);
	H5Pclose(plist_file_id);
	
	return 0;
}

int WriterParallel::overwrite_file_parallel(const std::string &filename)
{
	hid_t plist_file_id;
	MPI::Info info;

	// ==================================
	// Set up file access property list
	// ==================================
	plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_file_id, *comm_id_ptr, info);

	// ==================================
	// Create a new file collectively
	// ==================================
	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file_id);
	H5Pclose(plist_file_id);
	
	std::cout << "File created: " << file_id << std::endl;
	
	return 0;
}

int WriterParallel::writedata(long step, std::string const &dataset, const Parts &parts)
{
	// ==================================
	// Initialize all variables
	// ==================================
	MPI::Info info;

	long n_pts             = parts.n_pts;
	const double_vec * _x  = &parts.x;
	const double_vec * _xp = &parts.xp;
	const double_vec * _y  = &parts.y;
	const double_vec * _yp = &parts.yp;
	const double_vec * _z  = &parts.z;
	const double_vec * _zp = &parts.zp;
	double *buf;

	int n_write = ionsim::MAX_N_WRITE;

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
	// Access or create a new group
	// ==================================
	step_group_id    = group_step_access(file_id, step);

	// ==================================
	// Create dataset collectively
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	int rank = 2;
	count[0] = n_pts*p;
	count[1] = 6;
	dataset_id = dataset_create(step_group_id, dataspace_id, rank, count, dataset);

	// ==================================
	// Write dataset rows
	// ==================================
	plist_dx_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_dx_id, H5FD_MPIO_COLLECTIVE);

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
