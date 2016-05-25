#include <config.h>
#include <stdio.h>
#include <hdf5.h>
#include <string>
#include <iomanip>
#include "ions.h"
#include "support_func.h"
#include "emit.h"

void read_ion(hid_t step_group_id, std::string dataset)
{
	herr_t status;
	hid_t electrons_dataset_id;

	// ==================================
	// Data to read
	// ==================================
	int n_rows = 1;
	int n_cols = 6;
	int n_eles = n_rows*n_cols;
	double *buf;
	buf = new double[n_eles];

	// ==================================
	// Dataspace to read
	// ==================================
	hid_t electrons_dataspace_id;
	hsize_t start[2];
	start[0] = 0;
	start[1] = 0;
	hsize_t count[2];
	count[0] = n_rows;
	count[1] = n_cols;

	// ==================================
	// Memspace size
	// ==================================
	hsize_t buf_size[1];
	buf_size[0] = n_eles;

	// ==================================
	// Memspace to write to
	// ==================================
	hid_t memspace;
	hsize_t start_out[1];
	start_out[0] = 0;
	hsize_t count_out [1];
	count_out [0] = n_eles;

	electrons_dataset_id   = H5Dopen(step_group_id, dataset.c_str(), H5P_DEFAULT);

	// ==================================
	// Load electrons dataspace
	// ==================================
	electrons_dataspace_id = H5Dget_space(electrons_dataset_id);
	/* electrons_dataspace_id = H5Screate_simple(2, count, NULL); */
	status = H5Sselect_hyperslab(electrons_dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

	// ==================================
	// Load memspace dataspace
	// ==================================
	memspace = H5Screate_simple(1, buf_size, NULL);

	// ==================================
	// Attempt to read
	// ==================================
	status = H5Dread(electrons_dataset_id, H5T_NATIVE_DOUBLE, memspace, electrons_dataspace_id, H5P_DEFAULT, buf);

	std::cout << dataset << ": ";
	for (int i=0; i < n_eles; i++)
	{
		/* std::cout << buf[i]; */
		printf("% e", buf[i]);
		if (i+1 < n_eles)
		{
			std::cout << ", ";
		}
	}
	std::cout << std::endl;

	delete[] buf;
}

int main(int argc, char **argv)
{
	const char *filename;

	// ==================================
	// Determine whether to open default
	// ==================================
	std::cout << "=========================================" << std::endl;
	if (argc == 2)
	{
		filename = argv[1];
		std::cout << std::left << std::setw(25) << "Checking file: " << filename << std::endl;
	} else {
		filename = "output.h5";
		std::cout << std::left << std::setw(25) << "Checking default file: " << filename << std::endl;
	}
	std::cout << "=========================================" << std::endl;

	// ==================================
	// Variables for opening group
	// ==================================
	hid_t file_id;
	hid_t step_group_id;

	// ==================================
	// Load electrons dataset
	// ==================================
	file_id       = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	step_group_id = H5Oopen(file_id, "Step_0000/ions_steps", H5P_DEFAULT);

	read_ion(step_group_id, "ions_0000");
	read_ion(step_group_id, "ions_0001");

	// ==============================
	// Generate beam
	// ==============================
	long n_e                = 1e6;
	long n_ions             = 1e4;
	double q_tot            = 2e10;
	double radius           = 2.4276628847185805e-06;
	double sz               = 30e-6;
	double length           = 100e-6;
	double E                = 20.35;
	double emit_n           = 50e-6;
	double n_p_cgs          = 1e17;
	double m_ion_amu        = 1.00794;
	double sdelta           = 0.04;
	zdist_t zdist           = Z_DIST_FLAT;
	int n_steps             = 1;
	/* std::string filename    = "output.h5"; */
	pushmethod_t pushmethod = PUSH_SIMPLE;
	/* pushmethod_t pushmethod = PUSH_FIELD; */

	if (pushmethod == PUSH_SIMPLE)
	{
		filename = "simple.h5";
	} else if (pushmethod == PUSH_FIELD) {
		filename = "field.h5";
	} else {
		filename = "output.h5";
	}

	long n_field_x          = 101;
	long n_field_y          = 101;
	long n_field_z          = 51;
	double field_trans_wind = radius;

	double sr = ionsim::sr(emit_n, E, n_p_cgs, m_ion_amu);
	std::cout << "Sr: " << sr << std::endl;
	double nb_0 = ionsim::nb_0(q_tot, sz, sr);

	double z_end = (11.1367*GSL_CONST_MKSA_SPEED_OF_LIGHT / GSL_CONST_MKSA_ELECTRON_CHARGE) * sqrt(GSL_CONST_MKSA_VACUUM_PERMITTIVITY * m_ion_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS / nb_0);

	q_tot *= z_end / sz;
	sz = z_end;

	const SimParams simparams(
		E,
		emit_n,
		length,
		m_ion_amu,
		n_p_cgs,
		q_tot,
		radius,
		sz,
		sdelta,
		zdist,
		n_steps,
		pushmethod,
		n_e,
		n_field_x,
		n_field_y,
		n_field_z,
		field_trans_wind,
		z_end,
		n_ions,
		filename
		);
	std::cout << "F_r: " << F_r(2e-6, 0.0, simparams) << std::endl;

	// ==================================
	// Close everything
	// ==================================
	H5close();

	return 0;
}
