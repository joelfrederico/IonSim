#include <iomanip>
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
Ebeam::Ebeam(const SimParams &simparams, Beam x_beam, Beam y_beam, unsigned long int s) : Parts(simparams, PARTS_E), qpp(simparams.q_tot/n_pts), _simparams(simparams)
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

Ebeam::Ebeam(const SimParams &simparams, double sx, double sy, unsigned long int s) : Parts(simparams, PARTS_E), qpp(simparams.q_tot/n_pts), _simparams(simparams)
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
	double* rho_x;
	double* rho_y;
	double z_len;
	// ==================================
	// Set random number generator
	// ==================================
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	
	// ==================================
	// Set random number generator seed
	// ==================================
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
		
		z_len = sqrt(z_cov[0][0]);

		z[i] = gsl_ran_flat(r, -z_len/2, z_len/2);
		zp[i] = gsl_ran_gaussian(r, z_len);
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
	// Check limits are correct
	// ==================================
	if (z0 > z1)
	{
		std::swap(z1, z0);
	}

	// ==================================
	// Check and xfer particles
	// ==================================
	int reserve_length = 0;
	for (long i=0; i < n_pts; i++)
	{
		if ( (z0 < z[i]) && (z[i] < z1) ) {
			reserve_length++;
		}
	}

	x_out.reserve(reserve_length);
	xp_out.reserve(reserve_length);
	y_out.reserve(reserve_length);
	yp_out.reserve(reserve_length);
	z_out.reserve(reserve_length);
	zp_out.reserve(reserve_length);

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
			ind++;
		}
	}

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

int Ebeam::field_BE(Field_Data &field)
{
	double sx, sy, var_x, var_y, var_x_minus_var_y, f, temp, rsq, r;
	double xx, yy, x_in, y_in, xt, yt, Ex, Ey, Et, Em;
	bool sx_bigger;

	complex_double E;

	sx = x_std();
	sy = y_std();

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

			field.Ex_ind(i, j, k) = Ex;
			field.Ey_ind(i, j, k) = Ey;
		}
	}

	return 0;
}

int Ebeam::field_Coulomb(Field_Data &field)
{
	double dx, dy, dz;
	double drsq, dr, dr52;
	double x_e, y_e, z_e;

	std::cout << "Calculating Coulomb field..." << std::endl;

	const double common_para = qpp*GSL_CONST_MKSA_ELECTRON_CHARGE / (4*M_PI*GSL_CONST_MKSA_VACUUM_PERMITTIVITY);
	const double common_tran = common_para * _simparams.gamma_rel;
	double temp_para, temp_tran;

	for (int n=0; n < n_pts; n++) {
		x_e = x[n];
		y_e = y[n];
		z_e = z[n];

		for (int i=0; i < field.x_pts; i++) {
			dx = field.x_grid[i] - x_e;
			for (int j=0; j < field.y_pts; j++) {
				dy = field.y_grid[j] - y_e;
				for (int k=0; k < field.z_pts; k++) {
					dz = field.z_grid[k] - z_e;

					drsq = dx*dx + dy*dy + dz*dz;
					dr = sqrt(drsq);
					dr52 = drsq*dr;

					temp_para = common_para / dr52;
					temp_tran = common_tran / dr52;

					field.Ex_ind(i, j, k) = temp_tran * dx;
					field.Ey_ind(i, j, k) = temp_tran * dy;
					field.Ez_ind(i, j, k) = temp_para * dz;
				}
			}
		}
	}

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if (id == 1)
	{
		std::cout << "=============================================" << std::endl;

		std::cout << "dx: "           << dx                           << std::endl;
		std::cout << "dz: "           << dz                           << std::endl;
		std::cout << "z_e: "          << z_e                          << std::endl;
		std::cout << "z_grid: "       << field.z_grid[0]              << std::endl;
		std::cout << "gamma_rel: "    << _simparams.gamma_rel         << std::endl;
		std::cout << "eV: "           << GSL_CONST_MKSA_ELECTRON_VOLT << std::endl;
		std::cout << "E: "            << _simparams.E                 << std::endl;
		std::cout << ": "             << ELECTRON_REST_ENERGY         << std::endl;
		std::cout << "Stored field: " << std::setw(10)                << field.Ex_ind(4, 4, 1) << std::endl;

		std::cout << "=============================================" << std::endl;
	}

	return 0;
}
