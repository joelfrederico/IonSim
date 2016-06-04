#include "ebeam.h"
#include "hdf5_classes.h"
#include "loop_comm.h"
#include "support_func.h"
#include "writer_parallel.h"
#include <hdf5.h>
#include <mpi.h>

WriterParallel::WriterParallel(const std::string &filename, const MPI_Comm comm_id) :
	WriterBase(filename)
{
	_init(filename, false);
}

WriterParallel::WriterParallel(const std::string &filename, const MPI_Comm comm_id, bool overwrite) :
	WriterBase(filename)
{
	_init(filename, overwrite);
}

int WriterParallel::_init(const std::string &filename, bool overwrite)
{
	writer_type = WRITER_PARALLEL;
	if (overwrite) {
		overwrite_file_parallel();
	} else {
		open_file();
	}

	return 0;
}

int WriterParallel::open_file()
{
	// ==================================
	// Set up file access property list
	// ==================================
	MPI_Info info = MPI_INFO_NULL;

	PlistCreate plist_file(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_file.plist_id, loopcomm.slave_comm, info);

	// ==================================
	// Open the file
	// ==================================
	file_id = H5Fopen(_filename.c_str(), H5F_ACC_RDWR, plist_file.plist_id);
	
	return 0;
}

int WriterParallel::overwrite_file_parallel()
{
	MPI_Info info = MPI_INFO_NULL;

	// ==================================
	// Set up file access property list
	// ==================================
	PlistCreate plist_file(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_file.plist_id, loopcomm.slave_comm, info);

	// ==================================
	// Create a new file collectively
	// ==================================
	file_id = H5Fcreate(_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file.plist_id);
	
	return 0;
}

int WriterParallel::_writedata(DatasetAccess *&dataset, hid_t &loc_id, const std::string &dataset_str, const Parts &parts, unsigned int step)
{
	long long n_pts  = parts.n_pts;
	int n_write = MAX_N_WRITE;
	double *buf;

	const double_vec * _x  = &parts.x;
	const double_vec * _xp = &parts.xp;
	const double_vec * _y  = &parts.y;
	const double_vec * _yp = &parts.yp;
	const double_vec * _z  = &parts.z;
	const double_vec * _zp = &parts.zp;

	hsize_t count[2];
	hsize_t offset[2];
	hsize_t cdims[2];
	DataspaceCreate *memspace;
	herr_t status, status1;

	// ==================================
	// Create dataset collectively
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	int rank = 2;
	count[0] = n_pts*(loopcomm.p-1);
	count[1] = 6;

	// Creates new dataset_id, dataspace_id
	/* DatasetAccess dataset = DatasetAccess(loc_id, dataset_str, rank, count); */
	dataset = new DatasetAccess(loc_id, dataset_str, rank, count);

	AttributeCreate attribute(dataset->dataset_id, "Step", step);

	// ==================================
	// Write dataset rows
	// ==================================
	PlistCreate plist_dx = PlistCreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(plist_dx.plist_id, H5FD_MPIO_COLLECTIVE);

	PlistCreate dcpl(H5P_DATASET_CREATE);
	cdims[0] = 1000;
	cdims[1] = 6;
	H5Pset_chunk(dcpl.plist_id, 2, cdims);
	H5Pset_deflate(dcpl.plist_id, 6);

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

		offset[0] = (loopcomm.id-1) * n_pts + i;
		offset[1] = 0;


		status1 = H5Sselect_hyperslab(dataset->dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

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

		status = H5Dwrite(dataset->dataset_id, H5T_NATIVE_DOUBLE, memspace->dataspace_id, dataset->dataspace_id, plist_dx.plist_id, buf);
		if ((status < 0) && (loopcomm.id == 1))
		{
			std::cout << "==================================" << std::endl;

			std::cout << "PROBLEM"       << std::endl;
			std::cout << "status: "      << status                 << std::endl;
			std::cout << "status1: "     << status1                << std::endl;
			std::cout << "memspace_id: " << memspace->dataspace_id << std::endl;
			std::cout << "dataset_id: "  << dataset->dataset_id    << std::endl;
			std::cout << "dataset_str: " << dataset->_dataset_str  << std::endl;
			std::cout << "x[0]: "        << (*_x)[0]               << std::endl;
			std::cout << "buf[0]: "      << (*_x)[0]               << std::endl;

			std::cout << "==================================" << std::endl;
		}
		delete memspace;
	}
	delete [] buf;

	/* delete dataset; */

	return 0;
}

int WriterParallel::writedata(unsigned int step, const std::string &dataset_str, const Parts &parts)
{
	DatasetAccess *dataset;
	GroupStepAccess step_group = GroupStepAccess(file_id, step);

	_writedata(dataset, step_group.group_id, dataset_str, parts, step);

	delete dataset;

	return 0;
}

int WriterParallel::writedata(unsigned int step, const std::string &dataset_str, const Ebeam &ebeam)
{
	DatasetAccess *dataset;
	GroupStepAccess step_group = GroupStepAccess(file_id, step);

	_writedata(dataset, step_group.group_id, dataset_str, ebeam, step);

	AttributeCreate n_resolve ( dataset->dataset_id , "n_resolve" , ebeam.n_resolve() );
	AttributeCreate sr_macro  ( dataset->dataset_id , "sr_macro"  , ebeam.sr_macro()  );

	delete dataset;

	return 0;
}

int WriterParallel::writedata_substep(unsigned int step, unsigned int substep, const std::string &dataset_str, const std::string &subgroup, const Parts &parts)
{
	DatasetAccess *dataset;

	GroupStepAccess step_group = GroupStepAccess(file_id, step);
	GroupAccess sub_group = GroupAccess(step_group.group_id, subgroup);

	_writedata(dataset, sub_group.group_id, dataset_str, parts, substep);

	delete dataset;
	return 0;
}
