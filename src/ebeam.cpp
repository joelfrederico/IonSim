#include "ebeam.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "support_func.h"
#include <sstream>

// ==================================
// EBeam
// ==================================
Ebeam::Ebeam(int n_pts, double mass, double q_tot, double E, Beam x_beam, Beam y_beam, double z_cov[2][2]) : Parts(n_pts, mass)
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

int Ebeam::dump(std::string const &filename, int step, MPI::Intracomm &comm)
{
	std::stringstream dataset;
	dataset << step;
	ionsim::dump(filename, "ebeam", dataset.str(), comm, this);
	return 0;
}

