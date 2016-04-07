#include "ions.h"
#include "field_data.h"
#include "field_interp.h"
#include <complex>
#include <sstream>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include "support_func.h"
#include "simparams.h"
#include "loop_comm.h"

// ==============================
// Ions
// ==============================
Ions::Ions(const SimParams &simparams, Plasma &plasma, int n_pts, double radius, double length) : Parts(simparams, PARTS_ION)
{
	_plasma      = &plasma;
	_radius      = radius;
	/* printf("Received radius is: %0.6e\n", _radius); */
	/* printf("Received mass is: %0.6e\n", simparams.ion_mass()); */

	_part_charge = plasma.n_p() * length * M_PI * pow(radius, 2) * GSL_CONST_MKSA_ELECTRON_CHARGE;

	bool keep_looking;

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for (int i=0; i < n_pts; i++)
	{
		keep_looking = true;
		while (keep_looking)
		{
			x[i]     = gsl_ran_flat(r, -radius, radius);
			y[i]     = gsl_ran_flat(r, -radius, radius);
			if ((pow(x[i], 2.0) + pow(y[i], 2.0)) < pow(radius, 2.0))
			{
				keep_looking = false;
				if (x[i] > radius) printf("Radius: %0.6e", x[i]);
			}
		}
		z[i]  = gsl_ran_flat(r, 0, length);
		xp[i] = 0;
		yp[i] = 0;
		zp[i] = 0;
	}
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

	nb_0  = ((double *)params)[0];
	sig_r = ((double *)params)[1];
	mass  = ((double *)params)[2];

	x0  = &y[0];
	xp0 = &y[2];
	y0  = &y[1];
	yp0 = &y[3];


	F = F_r(*x0, *y0, nb_0, sig_r);

	dydt[0] = y[2];
	dydt[2] = F.real() / mass;
	dydt[1] = y[3];
	dydt[3] = F.imag() / mass;

	return GSL_SUCCESS;
}

int Ions::push(double dt, double nb_0, double sig_r)
{
	std::complex<double> F;
	double params[3] = {nb_0, sig_r, mass};
	double t = 0;

	gsl_odeiv2_system sys = {func, NULL, 4, &params};
	gsl_odeiv2_evolve *e  = gsl_odeiv2_evolve_alloc(4);
	gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, 4);

	/* gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, dt, sig_r/100, 1e-8); */
	
	for (int i=0; i < n_pts; i++)
	{
		double y[] = {x[i], y[i], xp[i], yp[i]};
		double yerr[4];
		/* gsl_odeiv2_driver_apply_fixed_step(d, &t, dt, 4, y); */
		gsl_odeiv2_step_apply(step, t, dt, y, yerr, NULL, NULL, &sys);
		x[i]  = y[0];
		xp[i] = y[2];
		y[i]  = y[1];
		yp[i] = y[3];

		F = F_r(x[i], y[i], nb_0, sig_r);
		z[i]  = std::abs(F);
		zp[i] = mass;
	}

	/* gsl_odeiv2_driver_free(d); */
	gsl_odeiv2_step_free(step);
	gsl_odeiv2_evolve_free(e);

	return 0;
}

int Ions::push_simple(double dt, double nb_0, double sig_r)
{
	std::complex<double> F;
	for (int i=0; i < n_pts; i++)
	{
		x[i] += xp[i] * dt;
		y[i] += yp[i] * dt;

		F = F_r(x[i], y[i], nb_0, sig_r);
		xp[i] += F.real() * dt / mass;
		yp[i] += F.imag() * dt / mass;

		z[i]  = F.real();
		zp[i] = F.imag();
	}
	return 0;
}

int Ions::push_field(double dt, Field_Data &field, int z_step)
{
	Field_Interp fieldinterp(field, *gsl_interp2d_bicubic);
	double Fx, Fy;

	for (int i=0; i < n_pts; i++)
	{
		Fx = fieldinterp.Ex(x[i], y[i], z_step);
		Fy = fieldinterp.Ey(x[i], y[i], z_step);

		x[i] += xp[i] * dt;
		y[i] += yp[i] * dt;

		xp[i] += Fx * dt / mass;
		yp[i] += Fy * dt / mass;

		z[i]  = Fx;
		zp[i] = Fy;
	}
	/* printf("Updated\n"); */
	return 0;
}
