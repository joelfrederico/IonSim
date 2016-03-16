#include "ebeam.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include "support_func.h"
#include <sstream>
#include "fields.h"
#include "consts.h"

// ==================================
// Constructors
// ==================================
Ebeam::Ebeam(SimParams &simparams, Beam x_beam, Beam y_beam) : Parts(simparams, ionsim::PARTS_E)
{
	_simparams = &simparams;
	q_tot      = simparams.q_tot;

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

	for (int i=0; i < n_pts(); i++)
	{
		rho_x = &x_cov[0][1];
		rho_y = &y_cov[0][1];

		gsl_ran_bivariate_gaussian(r , sqrt(x_cov[0][0]) , sqrt(x_cov[1][1]) , *rho_x , &x[i] , &xp[i]);
		gsl_ran_bivariate_gaussian(r , sqrt(y_cov[0][0]) , sqrt(y_cov[1][1]) , *rho_y , &y[i] , &yp[i]);
		/* gsl_ran_bivariate_gaussian(r , sqrt(z_cov[0][0]) , sqrt(z_cov[1][1]) , 0      , &_z[i] , &_zp[i]); */
		z[i] = gsl_ran_flat(r, 0, sqrt(z_cov[0][0]));
		zp[i] = gsl_ran_gaussian(r, sqrt(z_cov[1][1]));
	}
	gsl_rng_free(r);
}

Ebeam::Ebeam(
		SimParams &simparams,
		double_vec x_in,
 		double_vec xp_in,
 		double_vec y_in,
 		double_vec yp_in,
 		double_vec z_in,
 		double_vec zp_in
		) : Parts(simparams, ionsim::PARTS_E)
{
	_simparams = &simparams;

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
// Private Methods
// ==================================

double i_to_x(long i)
{
	return 0;
	/* long midpoint = */ 
}

double j_to_y(long j)
{
}

double k_to_z(long k)
{
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
	return gsl_stats_mean(x, 1, n_pts());
}

double Ebeam::x_std()
{
	double* x = &x[0];
	return gsl_stats_sd(x, 1, n_pts());
}

double Ebeam::y_mean()
{
	double* y = &y[0];
	return gsl_stats_mean(y, 1, n_pts());
}

double Ebeam::y_std()
{
	double* y = &y[0];
	return gsl_stats_sd(y, 1, n_pts());
}

Ebeam Ebeam::between(double z0, double z1)
{
	double_vec x_out;
	double_vec xp_out;
	double_vec y_out;
	double_vec yp_out;
	double_vec z_out;
	double_vec zp_out;

	long reserve_length = x.size();

	x_out.reserve(reserve_length);
	xp_out.reserve(reserve_length);
	y_out.reserve(reserve_length);
	yp_out.reserve(reserve_length);
	z_out.reserve(reserve_length);
	zp_out.reserve(reserve_length);

	if (z0 > z1)
	{
		std::swap(z1, z0);
	}

	long _n_pts = n_pts();
	long ind = 0;
	for (long i=0; i < _n_pts; i++)
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

	x.shrink_to_fit();
	xp.shrink_to_fit();
	y.shrink_to_fit();
	yp.shrink_to_fit();
	z.shrink_to_fit();
	zp.shrink_to_fit();

	return Ebeam(
			*(*this)._simparams,
			x,
			xp,
			y,
			yp,
			z,
			zp
		    );
}

int Ebeam::get_field(Field &field)
{
	complex_double E;
	double sx = x_std();
	double sy = y_std();
	double a, b;
	
	for (int i=0; i < field.x_pts; i++)
	{
		for (int j=0; j < field.y_pts; j++)
		{
			/* E = q_tot / (2* */
		}
	}
	return 0;
}
