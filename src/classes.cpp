#include "baseclass.h"
#include "classes.h"
#include "consts.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include "hdf5.h"

// ==============================
// Beam
// ==============================
Beam::Beam() {}

Beam::Beam(double beta, double alpha, Emit emit)
{
	_beta  = beta;
	_alpha = alpha;
	_emit = emit;
}

double Beam::alpha()
{
	return _alpha;
}

double Beam::beta()
{
	return _beta;
}

void Beam::cov(double output[2][2])
{
	output[0][0] = beta();
	output[0][1] = output[1][0] = -alpha();
	output[1][1] = (1.0+pow(alpha(), 2))/beta();
	double buf = _emit.emit();
	for (int i=0; i < 2; i++) 
	{
		for (int j=0; j < 2; j++)
		{
			output[i][j] *= _emit.emit();
		}
	}
	return;
}



// ==============================
// Ions
// ==============================
Ions::Ions() : Parts(0)
{
}

Ions::Ions(Plasma plasma, int n_pts, double x_window, double y_window) : Parts(n_pts)
{
	_plasma   = plasma;
	_x_window = x_window;
	_y_window = y_window;
	_n_pts    = n_pts;


	_delta.reserve(n_pts);

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	/* gsl_rng_set(r, MPI::COMM_WORLD.Get_rank() + 1); */

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for (int i=0; i < n_pts; i++)
	{
		_x[i]     = gsl_rng_uniform(r);
		_xp[i]    = gsl_rng_uniform(r);
		_y[i]     = gsl_rng_uniform(r);
		_yp[i]    = gsl_rng_uniform(r);
		_z[i]     = gsl_rng_uniform(r);
		_delta[i] = gsl_rng_uniform(r);
	}
	/* gsl_rng_free(r); */
}

double_vec Ions::x()
{
	return _x;
}

int Ions::dump(std::string const &filename)
{
	/* printf("Will print to: %s\n", filename.c_str()); */
	return 0;
}

// ==================================
// EBeam
// ==================================
Ebeam::~Ebeam() {}

Ebeam::Ebeam(int n_pts, double q_tot, double E, Beam x_beam, Beam y_beam, double z_cov[2][2]) : Parts(n_pts)
{
	_q_tot = q_tot;
	_E     = E;
	_x_beam = x_beam;
	_y_beam = y_beam;
	double* rho_x;
	double* rho_y;
	double x_cov[2][2];
	double y_cov[2][2];

	x_beam.cov(x_cov);
	y_beam.cov(y_cov);
	/* printf("x_cov: %.6d\n", x_cov); */

	gsl_rng * r = gsl_rng_alloc(gsl_rng_default);

	/* gsl_rng_env_setup(); */
	gsl_rng_set(r, MPI::COMM_WORLD.Get_rank() + 1);

	for (int i=0; i < n_pts; i++)
	{
		rho_x = &x_cov[0][1];
		rho_y = &y_cov[0][1];
		gsl_ran_bivariate_gaussian(r , sqrt(x_cov[0][0]) , sqrt(x_cov[1][1]) , *rho_x , &_x[i] , &_xp[i]);
		gsl_ran_bivariate_gaussian(r , sqrt(y_cov[0][0]) , sqrt(y_cov[1][1]) , *rho_y , &_y[i] , &_yp[i]);
		/* gsl_ran_bivariate_gaussian(r , sqrt(z_cov[0][0]) , sqrt(z_cov[1][1]) , 0      , &_z[i] , &_zp[i]); */
		_z[i] = gsl_ran_flat(r, 0, sqrt(z_cov[0][0]));
		_zp[i] = gsl_ran_gaussian(r, sqrt(z_cov[1][1]));
	}
	gsl_rng_free(r);
}

int Ebeam::dump(std::string const &filename, MPI::Intracomm &comm)
{
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
	memcpy(count, (hsize_t [2]){_n_pts*comm.Get_size(), 6}, 2*sizeof(hsize_t));
	dataspace_id = H5Screate_simple(2, count, NULL);

	// ==================================
	// Create dataset
	// ==================================
	dataset_id = H5Dcreate(group_id, "particles", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// ==================================
	// Write to hyperslab
	// ==================================
	memcpy(count, (hsize_t [2]){_n_pts, 6}, 2*sizeof(hsize_t));
	memcpy(offset, (hsize_t [2]){id * count[0], 0}, 2*sizeof(hsize_t));
	memspace_id = H5Screate_simple(2, count, NULL);

	H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

	plist_dx_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_dx_id, H5FD_MPIO_COLLECTIVE);
	double buf[_n_pts*6];
	for (int i=0; i < _n_pts; i++)
	{
		buf[i*6]     = _x[i];
		buf[i*6 + 1] = _xp[i];
		buf[i*6 + 2] = _y[i];
		buf[i*6 + 3] = _yp[i];
		buf[i*6 + 4] = _z[i];
		buf[i*6 + 5] = _zp[i];
	}

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, plist_dx_id, buf);

	// ==================================
	// Close file
	// ==================================
	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);
	H5Gclose(group_id);
	H5Fclose(file_id);

	return 0;
}

// ==================================
// Plasma
// ==================================
Plasma::Plasma() {}

Plasma::Plasma(double n_p_cgs, double ion_mass_amu)
{
	_n_p          = n_p_cgs * 1e6;
	_ion_mass_amu = ion_mass_amu;
}

double Plasma::n_p()
{
	return _n_p;
}

double Plasma::m()
{
	return _ion_mass_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
}

double Plasma::w_p()
{
	return sqrt( n_p() * pow(GSL_CONST_MKSA_ELECTRON_CHARGE, 2.0) / (GSL_CONST_MKSA_MASS_ELECTRON * GSL_CONST_MKSA_VACUUM_PERMITTIVITY) );
}

double Plasma::k_ion(double E)
{
	return n_p() * pow(GSL_CONST_MKSA_ELECTRON_CHARGE, 2.0) / (2.0 * E * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT * GSL_CONST_MKSA_VACUUM_PERMITTIVITY);
}

// ==============================
// Match
// ==============================
Match::Match(Plasma plasma, double E, Emit emit)
{
	_plasma = plasma;
	_E = E;
	_emit = emit;
}

double Match::beta()
{
	return 1.0 / sqrt( _plasma.k_ion(_E));
}

double Match::alpha()
{
	return 0;
}
