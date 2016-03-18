#include "ebeam.h"
#include "beam.h"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_exp.h>
#include "support_func.h"
#include <sstream>
#include "fields.h"
#include "consts.h"
#include <math.h>
#include "simparams.h"
#include "faddeeva/Faddeeva.hh"

// ==================================
// Constructors
// ==================================
Ebeam::Ebeam(SimParams &simparams, Beam x_beam, Beam y_beam) : Parts(simparams, ionsim::PARTS_E), qpp(simparams.q_tot/n_pts), _simparams(simparams)
{
	// ==================================
	// Save things
	// ==================================
	_x_beam = x_beam;
	_y_beam = y_beam;

	// ==================================
	// Initialize local variables
	// ==================================
	double* rho_x;
	double* rho_y;
	double x_cov[2][2];
	double y_cov[2][2];

	// ==================================
	// Get covariances for x and y
	// ==================================
	x_beam.cov(x_cov);
	y_beam.cov(y_cov);

	// ==================================
	// Set random number generator
	// ==================================
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	
	// ==================================
	// Set random number generator seed
	// ==================================
	gsl_rng_set(r, MPI::COMM_WORLD.Get_rank() + 1);

	// ==================================
	// Create particles randomly
	// ==================================
	for (int i=0; i < n_pts; i++)
	{
		rho_x = &x_cov[0][1];
		rho_y = &y_cov[0][1];

		gsl_ran_bivariate_gaussian(r , sqrt(x_cov[0][0]) , sqrt(x_cov[1][1]) , *rho_x , &x[i] , &xp[i]);
		gsl_ran_bivariate_gaussian(r , sqrt(y_cov[0][0]) , sqrt(y_cov[1][1]) , *rho_y , &y[i] , &yp[i]);
		/* gsl_ran_bivariate_gaussian(r , sqrt(z_cov[0][0]) , sqrt(z_cov[1][1]) , 0      , &_z[i] , &_zp[i]); */
		z[i] = gsl_ran_flat(r, 0, sqrt(z_cov[0][0]));
		zp[i] = gsl_ran_gaussian(r, sqrt(z_cov[1][1]));
	}

	// ==================================
	// Free memory
	// ==================================
	gsl_rng_free(r);
}

Ebeam::Ebeam(
		SimParams &simparams,
		const double n_pts,
		const double type,
		double_vec x_in,
 		double_vec xp_in,
 		double_vec y_in,
 		double_vec yp_in,
 		double_vec z_in,
 		double_vec zp_in
		) : Parts(simparams.ion_mass(), n_pts, type), qpp(simparams.q_tot/n_pts), _simparams(simparams)
{
	x_in.shrink_to_fit();
	xp_in.shrink_to_fit();
	y_in.shrink_to_fit();
	yp_in.shrink_to_fit();
	z_in.shrink_to_fit();
	zp_in.shrink_to_fit();

	x  = x_in;
	xp = xp_in;
	y  = y_in;
	yp = yp_in;
	z  = z_in;
	zp = zp_in;
}

// ==================================
// Public Methods
// ==================================

int Ebeam::dump(std::string const &filename, int step, MPI::Intracomm &comm)
{
	std::stringstream dataset;
	dataset << step;
	ionsim::dump(filename, "ebeam", dataset.str(), comm, *this);
	return 0;
}

double Ebeam::x_mean()
{
	double* x = &x[0];
	return gsl_stats_mean(x, 1, n_pts);
}

double Ebeam::x_std()
{
	double* x = &x[0];
	return gsl_stats_sd(x, 1, n_pts);
}

double Ebeam::y_mean()
{
	double* y = &y[0];
	return gsl_stats_mean(y, 1, n_pts);
}

double Ebeam::y_std()
{
	double* y = &y[0];
	return gsl_stats_sd(y, 1, n_pts);
}

Ebeam Ebeam::between(double z0, double z1)
{
	// ==================================
	// Initialize vectors to hold coords
	// ==================================
	double_vec x_out;
	double_vec xp_out;
	double_vec y_out;
	double_vec yp_out;
	double_vec z_out;
	double_vec zp_out;

	// ==================================
	// Reserve space for all possible 
	// particles
	// ==================================
	long reserve_length = x.size();

	x_out.reserve(reserve_length);
	xp_out.reserve(reserve_length);
	y_out.reserve(reserve_length);
	yp_out.reserve(reserve_length);
	z_out.reserve(reserve_length);
	zp_out.reserve(reserve_length);

	// ==================================
	// Check limits are correct
	// ==================================
	if (z0 > z1)
	{
		std::swap(z1, z0);
	}

	// ==================================
	// Check and xfer particles
	// ==================================
	long ind = 0;
	for (long i=0; i < n_pts; i++)
	{
		if ( (z0 < z[i]) && (z[i] < z1) )
		{
			x_out[ind]  = x[i];
			xp_out[ind] = xp[i];
			y_out[ind]  = y[i];
			yp_out[ind] = yp[i];
			z_out[ind]  = z[i];
			zp_out[ind] = zp[i];
		}
	}

	// ==================================
	// Shrink arrays to reduce memory
	// ==================================
	x_out.shrink_to_fit();
	xp_out.shrink_to_fit();
	y_out.shrink_to_fit();
	yp_out.shrink_to_fit();
	z_out.shrink_to_fit();
	zp_out.shrink_to_fit();

	// ==================================
	// Create new ebeam from particles
	// ==================================
	return Ebeam(
			(*this)._simparams,
			x_out.size(),
			(*this).type,
			x_out,
			xp_out,
			y_out,
			yp_out,
			z_out,
			zp_out
		    );
}

int Ebeam::field(Field &field)
{
	complex_double E;
	double sx = x_std();
	double sy = y_std();

	bool sx_bigger = (sx > sy);
	if (!sx_bigger)
	{
		std::swap(sx, sy);
	}

	double var_x = pow(sx, 2);
	double var_y = pow(sx, 2);
	double a, b, f;

	double xx, yy, xt, yt, Ex, Ey, Et;
	double var_x_minus_var_y = var_x-var_y;
	
	for (long i=0; i < field.x_pts; i++)
	{
		xt = field.i_to_x(i);
		if (sx_bigger)
		{
			xx = xt;
		} else {
			yy = xt;
		}

		for (long j=0; j < field.y_pts; j++)
		{
			yt = field.j_to_y(j);
			if (sx_bigger)
			{
				yy = yt;
			} else {
				xx = -yt;
			}

			f = sqrt(2*var_x_minus_var_y);
			b = yy / f;
			a = xx / f;
			E = qpp*n_pts / (2*GSL_CONST_MKSA_VACUUM_PERMITTIVITY*sqrt(2*M_PI*var_x_minus_var_y)) * ( Faddeeva::w(complex_double(xx, yy)/f) - gsl_sf_exp(-pow(xx,2)/(2*var_x)+pow(yy,2)/(2*var_y))* Faddeeva::w(complex_double(xx*sy/sx, yy*sx/sy)/f) );
			Ex = E.imag();
			Ey = E.real();
			if (!sx_bigger) 
			{ 
				Et = Ex;
				Ex = Ey;
				Ey = -Et;
			}

			field.Ex_ind(i, j) = Ex;
			field.Ey_ind(i, j) = Ey;
		}
	}
	return 0;
}
