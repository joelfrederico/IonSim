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
#include "support_func.h"

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

Ions::Ions(Plasma * plasma, int n_pts, double radius, double length) : Parts(n_pts)
{
	_plasma      = plasma;
	_radius      = radius;
	_part_charge = (*plasma).n_p() * length * M_PI * pow(radius, 2) * GSL_CONST_MKSA_ELECTRON_CHARGE;

	bool keep_looking;

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	/* gsl_rng_set(r, MPI::COMM_WORLD.Get_rank() + 1); */

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for (int i=0; i < n_pts; i++)
	{
		keep_looking = true;
		while (keep_looking)
		{
			_x[i]     = gsl_ran_flat(r, -radius, radius);
			_y[i]     = gsl_ran_flat(r, -radius, radius);
			if ((pow(_x[i], 2) + pow(_y[i], 2)) < pow(radius, 2)) keep_looking = false;
		}
		_z[i]  = gsl_ran_flat(r, 0, length);
		_xp[i] = 0;
		_yp[i] = 0;
		_zp[i] = 0;
	}
}

int Ions::dump(std::string const &filename, MPI::Intracomm &comm)
{
	ionsim::dump(filename, comm, this);
	return 0;
}

// ==================================
// EBeam
// ==================================
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

	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

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
	ionsim::dump(filename, comm, this);
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
