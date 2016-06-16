#include "beam.h"
#include "consts.h"
#include "ebeam.h"
#include "faddeeva/Faddeeva.hh"
#include "field_data.h"
#include "gsl_classes.h"
#include "simparams.h"
#include "support_func.h"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_statistics_double.h>
#include <iomanip>
#include <math.h>
#include <mpi.h>
#include <gflags/gflags.h>

DECLARE_bool(verbose);

Cov z_cov(SimParams simparams)
{
	Cov out;
	out(0, 0) = simparams.sz * simparams.sz;
	out(0, 1) = out(1, 0) = 0;
	out(1, 1) = simparams.sdelta;

	return out;
}

// ==================================
// Constructors
// ==================================
Ebeam::Ebeam(const SimParams &simparams, Beam x_beam, Beam y_beam, unsigned long int s) : Parts(simparams, PARTS_E), _simparams(simparams), _x_cov(x_beam.cov()), _y_cov(y_beam.cov()), _z_cov(z_cov(simparams))
{
	_gen_bivariate_gaussian(s, _x_cov, _y_cov, _z_cov);
}

Ebeam::Ebeam(const SimParams &simparams, double sx, double sxp, double sy, double syp, unsigned long int s) : Parts(simparams, PARTS_E), _simparams(simparams), _x_cov(sx*sx, 0, 0, syp*syp), _y_cov(sy*sy, 0, 0, syp*syp), _z_cov(z_cov(simparams))
{
	_gen_bivariate_gaussian(s, _x_cov, _y_cov, _z_cov);
}

int Ebeam::_gen_bivariate_gaussian(unsigned long int s, Cov x_cov, Cov y_cov, Cov z_cov)
{
	int status;
	double z_temp;
	double* rho_x;
	double* rho_y;
	double z_len;
	int n_e_node;
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	// ==================================
	// Set random number generator
	// ==================================
	GSLrng rng(s);
	
	// ==================================
	// Create particles randomly
	// ==================================
	rho_x = &x_cov(0, 1);
	rho_y = &y_cov(0, 1);

	n_e_node = _simparams.n_e_node();

	if (FLAGS_verbose && (id == 1)) std::cout << "Creating n_pts: " << n_e_node << std::endl;

	for (int i=0; i < n_e_node; i++)
	{

		gsl_ran_bivariate_gaussian(rng.r , sqrt(x_cov(0, 0)) , sqrt(x_cov(1, 1)) , *rho_x , &x[i] , &xp[i]);
		gsl_ran_bivariate_gaussian(rng.r , sqrt(y_cov(0, 0)) , sqrt(y_cov(1, 1)) , *rho_y , &y[i] , &yp[i]);
		
		switch (_simparams.zdist)
		{
			case Z_DIST_FLAT:
				z_len = _simparams.sz;

				z[i] = gsl_ran_flat(rng.r, _simparams.z_center, _simparams.sz);

				break;

			case Z_DIST_GAUSS:
				gsl_ran_bivariate_gaussian(rng.r, sqrt(z_cov(0, 0)), sqrt(z_cov(1, 1)), 0, &z_temp, &zp[i]);
				z[i] = z_temp + _simparams.z_center;

				break;
		}

		zp[i] = gsl_ran_gaussian(rng.r, _simparams.sdelta);
	}

	// ==================================
	// Set resolution
	// ==================================

	if (verbose())
	{
		std::cout << "=========================================" << std::endl;
		std::cout << "Ebeam Calculated"                          << std::endl;
		std::cout << "-----------------------------------------" << std::endl;
		std::cout << "n_0: "       << n_0()                      << std::endl;
		std::cout << "n_resolve: " << n_resolve()                << std::endl;
		std::cout << "=========================================" << std::endl;
	}
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
		) : Parts(simparams.ion_mass(), simparams.n_e_node(), type), _simparams(simparams)
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
	for (long long i=0; i < n_pts; i++)
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

	long long ind = 0;
	for (long long i=0; i < n_pts; i++)
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
	temp = _simparams.qpp_e()*n_pts * GSL_CONST_MKSA_ELECTRON_CHARGE/GSL_CONST_MKSA_VACUUM_PERMITTIVITY;

	int k = 0;

	for (int i=0; i < field.x_pts; i++)
	{
		xt = field.x_grid[i];
		if (sx_bigger)
		{
			xx = xt;
		} else {
			yy = xt;
		}

		for (int j=0; j < field.y_pts; j++)
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

	const double common_para = _simparams.qpp_e()*GSL_CONST_MKSA_ELECTRON_CHARGE / (4*M_PI*GSL_CONST_MKSA_VACUUM_PERMITTIVITY);
	const double common_tran = common_para * _simparams.gamma_rel();
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

					field.Ex_ind(i, j, k) += temp_tran * dx;
					field.Ey_ind(i, j, k) += temp_tran * dy;
					field.Ez_ind(i, j, k) += temp_para * dz;
				}
			}
		}
	}

	return 0;
}

int Ebeam::field_Coulomb_sliced(Field_Data &field)
{
	double dx, dy;
	double drsq, dr;
	double x_e, y_e, z_e;
	const double temp = 0.5;
	double k_double;
	int k;
	double sr_m;
	double srsq_macro;
	const double dz = _simparams.dz();

	sr_m = sr_macro();
	srsq_macro = sr_m*sr_m;

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	/* std::cout << "Id: " << id << ", Qpp: " << _simparams.qpp_e() << std::endl; */
	const double common   = _simparams.qpp_e()*GSL_CONST_MKSA_ELECTRON_CHARGE / (4*M_PI*dz*srsq_macro);

	const double common_E = common / GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	const double common_B = common * GSL_CONST_MKSA_SPEED_OF_LIGHT * GSL_CONST_MKSA_VACUUM_PERMEABILITY;

	double temp_tran;
	double temp_tran_E;
	double temp_tran_B;

	for (int n=0; n < n_pts; n++) {
		x_e = x[n];
		y_e = y[n];
		z_e = z[n];

		k_double = z_e / dz;
		k = floor(k_double);
		field.Ez_ind(0, 1, 0)++;
		if (k < 0)
		{
			std::cout << "K: " << k << std::endl;
		}

		if ((0 <= k) && (k <= field.z_pts-1) )
		{
			field.Ez_ind(k, k, 0) ++;
			field.Ez_ind(0, 2, 0) ++;

			field.Bz_ind(k, k, 0) ++;
			/* field.Bz_ind(1, 3, 0) ++; */

			for (int i=0; i < field.x_pts; i++)
			{
				dx = field.x_grid[i] - x_e;
				for (int j=0; j < field.y_pts; j++)
				{
					dy = field.y_grid[j] - y_e;

					drsq = dx*dx + dy*dy;
					dr   = sqrt(drsq);

					temp_tran = gsl_sf_exprel(- drsq / (2*srsq_macro));

					temp_tran_E = temp_tran * common_E;
					temp_tran_B = temp_tran * common_B;

					field.Ex_ind(i, j, k) += temp_tran_E * dx;
					field.Ey_ind(i, j, k) += temp_tran_E * dy;

					field.Bx_ind(i, j, k) += -temp_tran_B * dy;
					field.By_ind(i, j, k) +=  temp_tran_B * dx;
				}
			}
		} else {
			field.Ez_ind(1, 0, 0)++;
			std::cout << "K: " << k << ", k_double: " << k_double << std::endl;
			std::cout << "zend: " << _simparams.z_end << ", dz: " << dz << ", z_e: " << z_e << std::endl;
		}
	}


	return 0;
}

int _calc_n_0_n_resolve(Cov x_cov, SimParams _simparams, double &n_0, double &n_resolve)
{
	int status;
	double srsq = x_cov(0, 0);

	switch (_simparams.zdist)
	{
		case Z_DIST_FLAT:
			n_0 = (_simparams.n_e) / (M_PI * srsq * _simparams.sz);
			n_resolve = n_0;
			/* n_resolve = -(_simparams.n_e) / (pow(2*M_PI, 1.5) * srsq * z_len * gsl_sf_expm1(-4.5)); */

			break;

		case Z_DIST_GAUSS:
			gsl_sf_result result;
			n_0 = (_simparams.n_e) / (pow(2*M_PI, 1.5) * srsq * _simparams.sz);
			/* n_resolve = gsl_sf_exp_mult(-4.5, n_0); */
			status = gsl_sf_exp_mult_e(-4.5, n_0, &result);

			if (status != GSL_SUCCESS)
			{
				std::cout << "Exp_mult_e error. n_0: " << n_0 << std::endl;
				return -1;
			}

			break;
	}

	return 0;
}

double Ebeam::n_0() const
{
	double n_0, n_resolve;
	_calc_n_0_n_resolve(_x_cov, _simparams, n_0, n_resolve);
	return n_0;
}

double Ebeam::n_resolve() const
{
	double n_0, n_resolve;
	_calc_n_0_n_resolve(_x_cov, _simparams, n_0, n_resolve);
	return n_resolve;
}

double Ebeam::sr_macro() const
{
	return 0.2347535410605456 / sqrt(n_0() * _simparams.dz());
}

bool Ebeam::verbose()
{
	int id;
	bool out;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if (FLAGS_verbose && (id == 1))
	{
		out = true;
	} else {
		out = false;
	}

	return out;
}
