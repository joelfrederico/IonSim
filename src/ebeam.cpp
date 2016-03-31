#include "ebeam.h"
#include <mpi.h>
#include "beam.h"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_exp.h>
#include "support_func.h"
#include <sstream>
#include "field_data.h"
#include "consts.h"
#include <math.h>
#include "simparams.h"
#include "faddeeva/Faddeeva.hh"

// ==================================
// Constructors
// ==================================
Ebeam::Ebeam(const SimParams &simparams, Beam x_beam, Beam y_beam, unsigned long int s) : Parts(simparams, ionsim::PARTS_E), qpp(simparams.q_tot/n_pts), _simparams(simparams)
{
	// ==================================
	// Save things
	// ==================================
	/* _x_beam = x_beam; */
	/* _y_beam = y_beam; */

	// ==================================
	// Initialize local variables
	// ==================================
	double x_cov[2][2];
	double y_cov[2][2];
	double z_cov[2][2];

	// ==================================
	// Get covariances for x and y
	// ==================================
	x_beam.cov(x_cov);
	y_beam.cov(y_cov);

	z_cov[0][0] = simparams.sz;
	z_cov[0][1] = z_cov[1][0] = 0;
	z_cov[1][1] = simparams.sdelta;

	_gen_bivariate_gaussian(s, x_cov, y_cov, z_cov);
}

Ebeam::Ebeam(const SimParams &simparams, double sx, double sy, unsigned long int s) : Parts(simparams, ionsim::PARTS_E), qpp(simparams.q_tot/n_pts), _simparams(simparams)
{
	// ==================================
	// Initialize local variables
	// ==================================
	double x_cov[2][2];
	double y_cov[2][2];
	double z_cov[2][2];

	// ==================================
	// Get covariances for x and y
	// ==================================
	for (int i=0; i < 2; i++)
	{
		for (int j=0; j < 2; j++)
		{
			x_cov[i][j] = 0;
			y_cov[i][j] = 0;
		}
	}
	x_cov[0][0] = sx;
	y_cov[0][0] = sy;

	z_cov[0][0] = simparams.sz;
	z_cov[0][1] = z_cov[1][0] = 0;
	z_cov[1][1] = simparams.sdelta;

	_gen_bivariate_gaussian(s, x_cov, y_cov, z_cov);
}

int Ebeam::_gen_bivariate_gaussian(unsigned long int s, double x_cov[2][2], double y_cov[2][2], double z_cov[2][2])
{
	/* std::cout << sqrt(x_cov[0][0]) << std::endl; */
	/* std::cout << sqrt(x_cov[1][1]) << std::endl; */
	/* std::cout << sqrt(y_cov[0][0]) << std::endl; */
	/* std::cout << sqrt(y_cov[1][1]) << std::endl; */
	
	double* rho_x;
	double* rho_y;
	// ==================================
	// Set random number generator
	// ==================================
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	
	// ==================================
	// Set random number generator seed
	// ==================================
	/* gsl_rng_set(r, MPI::COMM_WORLD.Get_rank() + 1); */
	gsl_rng_set(r, s);

	// ==================================
	// Create particles randomly
	// ==================================
	rho_x = &x_cov[0][1];
	rho_y = &y_cov[0][1];
	for (int i=0; i < n_pts; i++)
	{

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
	return 0;
}

Ebeam::Ebeam(
		const SimParams &simparams,
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

int Ebeam::dump_parallel(std::string const &filename, int step, MPI::Intracomm &comm)
{
	ionsim::dump_parallel(filename, step, "ebeam", comm, *this);
	return 0;
}

int Ebeam::dump_serial(std::string const &filename, int step)
{
	std::string group = "ebeam";
	std::string dataset = "particles";
	ionsim::dump_serial(filename, step, group, dataset, *this);
	return 0;
}

double Ebeam::x_mean()
{
	double* x = &x[0];
	return gsl_stats_mean(x, 1, n_pts);
}

double Ebeam::x_std()
{
	/* double* x_local = &x[0]; */
	return gsl_stats_sd_m(x.data(), 1, n_pts, 0);
}

double Ebeam::y_mean()
{
	double* y = &y[0];
	return gsl_stats_mean(y, 1, n_pts);
}

double Ebeam::y_std()
{
	/* double* y_local = &y[0]; */
	return gsl_stats_sd_m(y.data(), 1, n_pts, 0);
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

int Ebeam::field(Field_Data &field)
{
	double sx, sy, var_x, var_y, var_x_minus_var_y, f, temp, rsq, r;
	double xx, yy, x_in, y_in, xt, yt, Ex, Ey, Et, Em;
	bool sx_bigger;

	complex_double E;

	sx = x_std();
	sy = y_std();
	/* std::cout << "Sx: " << sx << ", Sy: " << sy << std::endl; */

	sx_bigger = (sx > sy);
	if (!sx_bigger)
	{
		std::swap(sx, sy);
	}

	var_x = sx*sx;
	var_y = sy*sy;

	var_x_minus_var_y = var_x-var_y;
	if (var_x_minus_var_y/var_x == 0) std::cout << "Using symmetric fields" << std::endl;
	temp = qpp*n_pts * GSL_CONST_MKSA_ELECTRON_CHARGE/GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	/* temp = qpp*n_pts / (2*GSL_CONST_MKSA_VACUUM_PERMITTIVITY*sqrt(2*M_PI*var_x_minus_var_y)); */
	/* int id = MPI::COMM_WORLD.Get_rank(); */
	/* printf("Proc: %d, Field with associated charge is: %0.3e\n", id, temp); */
	int k = 0;

	for (long i=0; i < field.x_pts; i++)
	{
		xt = field.x_grid[i];
		if (sx_bigger)
		{
			xx = xt;
		} else {
			yy = xt;
		}

		for (long j=0; j < field.y_pts; j++)
		{
			yt = field.y_grid[j];
			if (sx_bigger)
			{
				yy = yt;
			} else {
				xx = -yt;
			}

			rsq = xx*xx + yy*yy;
			if (rsq == 0) {
				Ex = 0;
				Ey = 0;
			} else if (var_x_minus_var_y/var_x < 1e-3) {
				r = sqrt(rsq);
				/* Et = -(qpp*n_pts / (GSL_CONST_MKSA_VACUUM_PERMITTIVITY)); */
				/* Et = 1; */
				/* Et = temp * -gsl_sf_exprel(-rsq/(2*var_x)) * r / (4*M_PI*var_x); */
				Et = temp * (1-exp(-rsq/(2*var_x)))/(2*M_PI*r);
				Ex = Et*xx/r;
				Ey = Et*yy/r;
			} else {
				y_in = std::abs(yy);
				f = sqrt(2*var_x_minus_var_y);
				E =  temp * ( Faddeeva::w(complex_double(xx, y_in)/f) - exp(-xx*xx/(2*var_x)-y_in*y_in/(2*var_y))* Faddeeva::w(complex_double(xx*sy/sx, y_in*sx/sy)/f) ) / (2* sqrt(2*M_PI * var_x_minus_var_y));

				Ex = E.imag();
				Ey = E.real()*GSL_SIGN(yy);

			}

			if (!sx_bigger) 
			{ 
				Et = Ex;
				Ex = Ey;
				Ey = -Et;
			}
			/* Ex = sqrt(Ex*Ex+Ey*Ey); */

			field.Ex_ind(i, j, k) = Ex;
			field.Ey_ind(i, j, k) = Ey;
		}
	}
	/* field.Ex_ind(25, 10) = 1; */
	return 0;
}
