#include "ions.h"
#include <complex>
#include <sstream>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include "support_func.h"

// ==============================
// Ions
// ==============================
Ions::Ions() : Parts(0, 0)
{
}

Ions::Ions(Plasma * plasma, int n_pts, double radius, double length) : Parts(n_pts, (*plasma).m())
{
	_plasma      = plasma;
	_radius      = radius;
	printf("Received radius is: %0.6e\n", _radius);
	printf("Received mass is: %0.6e\n", _mass);

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
			if ((pow(_x[i], 2.0) + pow(_y[i], 2.0)) < pow(radius, 2.0))
			{
				keep_looking = false;
				if (_x[i] > radius) printf("Radius: %0.6e", _x[i]);
			}
		}
		_z[i]  = gsl_ran_flat(r, 0, length);
		_xp[i] = 0;
		_yp[i] = 0;
		_zp[i] = 0;
	}
}

int Ions::dump(std::string const &filename, int step, MPI::Intracomm &comm)
{
	std::stringstream dataset;
	dataset << step;
	ionsim::dump(filename, "ions", dataset.str(), comm, this);
	return 0;
}

std::complex<double> F_r(double x, double y, double nb_0, double sig_r)
{
	double r = gsl_hypot(x, y);
	double Er_factor;
	double Er;
	double theta;

	Er_factor = GSL_CONST_MKSA_ELECTRON_CHARGE*GSL_CONST_MKSA_ELECTRON_CHARGE*nb_0*sig_r*sig_r / (2*M_PI*GSL_CONST_MKSA_VACUUM_PERMITTIVITY);
	Er = Er_factor * gsl_expm1(-r*r / (2*sig_r*sig_r)) / r;

	return std::complex<double>(Er * x/r, Er* y/r);
}

int func(double t, const double y[], double dydt[], void * params)
{
	(void)(t);
	double nb_0, sig_r, mass;
	std::complex<double> F;
	const double *x0, *xp0, *y0, *yp0;
	int p = MPI::COMM_WORLD.Get_rank();

	nb_0  = ((double *)params)[0];
	sig_r = ((double *)params)[1];
	mass  = ((double *)params)[2];
	/* printf("Herro: %d, %0.3e; %0.3e; %0.3e\n", p, nb_0, sig_r, mass); */

	x0  = &y[0];
	xp0 = &y[2];
	y0  = &y[1];
	yp0 = &y[3];


	F = F_r(*x0, *y0, nb_0, sig_r);

	dydt[0] = y[2];
	dydt[2] = F.real() / mass;
	dydt[1] = y[3];
	dydt[3] = F.imag() / mass;

	/* printf("%d alive!\n", p); */

	/* double *arr = *(double *)params; */
	return GSL_SUCCESS;
}

int Ions::push(double dt, double nb_0, double sig_r)
{
	std::complex<double> F;
	double params[3] = {nb_0, sig_r, _mass};
	double t = 0;
	int p = MPI::COMM_WORLD.Get_rank();

	gsl_odeiv2_system sys = {func, NULL, 4, &params};
	gsl_odeiv2_evolve *e  = gsl_odeiv2_evolve_alloc(4);
	gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, 4);

	/* gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, dt, sig_r/100, 1e-8); */


	for (int i=0; i < _n_pts; i++)
	{
		double y[] = {_x[i], _y[i], _xp[i], _yp[i]};
		double yerr[4];
		/* printf("%d Trying\n", p); */
		/* gsl_odeiv2_driver_apply_fixed_step(d, &t, dt, 4, y); */
		gsl_odeiv2_step_apply(step, t, dt, y, yerr, NULL, NULL, &sys);
		/* printf("%d Got here\n", p); */
		_x[i]  = y[0];
		_xp[i] = y[2];
		_y[i]  = y[1];
		_yp[i] = y[3];

		F = F_r(_x[i], _y[i], nb_0, sig_r);
		_z[i]  = std::abs(F);
		_zp[i] = _mass;
	}

	/* gsl_odeiv2_driver_free(d); */
	gsl_odeiv2_step_free(step);
	gsl_odeiv2_evolve_free(e);

	return 0;
}

int Ions::push_simple(double dt, double nb_0, double sig_r)
{
	std::complex<double> F;
	for (int i=0; i < _n_pts; i++)
	{
		_x[i] += _xp[i] * dt;
		_y[i] += _yp[i] * dt;

		F = F_r(_x[i], _y[i], nb_0, sig_r);
		_xp[i] += F.real() * dt / _mass;
		_yp[i] += F.imag() * dt / _mass;
	}
	return 0;
}

