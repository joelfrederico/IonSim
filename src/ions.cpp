#include "ions.h"
#include "field_data.h"
#include "field_interp.h"
#include <complex>
#include <sstream>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_exp.h>
#include "support_func.h"
#include "simparams.h"
#include "loop_comm.h"

// ==============================
// Ions
// ==============================
Ions::Ions(const SimParams *simparams, Plasma &plasma, int n_pts, double radius, double length) : Parts(*simparams, PARTS_ION), _simparams(simparams)
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

	x[0]  = 2e-6;
	xp[0] = 0;
	y[0]  = 0;
	yp[0] = 0;
	z[0]  = 0;
	zp[0] = 0;

	for (int i=1; i < n_pts; i++)
	/* for (int i=0; i < n_pts; i++) */
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
		/* z[i]  = gsl_ran_flat(r, 0, length); */
		z[i] = 0;
		xp[i] = 0;
		yp[i] = 0;
		zp[i] = 0;
	}
}

std::complex<double> F_r(double x, double y, const SimParams &simparams)
{
	double Er_factor;
	double Fr;
	double theta;

	double sr = ionsim::sr(simparams.emit_n, simparams.E, simparams.n_p_cgs, simparams.m_ion_amu);

	Er_factor = simparams.q_tot * GSL_CONST_MKSA_ELECTRON_CHARGE / (4*M_PI*GSL_CONST_MKSA_VACUUM_PERMITTIVITY*simparams.sz*sr*sr);
	Fr = -GSL_CONST_MKSA_ELECTRON_CHARGE * Er_factor * gsl_sf_exprel(-(x*x+y*y) / (2*sr*sr));

	/* Er = (simparams.q_tot * GSL_CONST_MKSA_ELECTRON_CHARGE) / (4 * M_PI * GSL_CONST_MKSA_VACUUM_PERMITTIVITY * simparams.sz * sr * sr); */

	return std::complex<double>(Fr * x, Fr* y);
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


	/* F = F_r(*x0, *y0, nb_0, sig_r); */

	dydt[0] = y[2];
	dydt[2] = F.real() / mass;
	dydt[1] = y[3];
	dydt[3] = F.imag() / mass;

	return GSL_SUCCESS;
}

int Ions::push(double nb_0, double sig_r)
{
	double dt = _simparams->dt();
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

		/* F = F_r(x[i], y[i], nb_0, sig_r); */
		z[i]  = std::abs(F);
		zp[i] = mass;
	}

	/* gsl_odeiv2_driver_free(d); */
	gsl_odeiv2_step_free(step);
	gsl_odeiv2_evolve_free(e);

	return 0;
}

int Ions::push_simple(double nb_0, double sig_r)
{
	double dt = _simparams->dt();
	std::complex<double> F;
	for (int i=0; i < n_pts; i++)
	{
		F = F_r(x[i], y[i], *_simparams);
		xp[i] = xp[i] + F.real() * dt / mass;
		yp[i] = yp[i] + F.imag() * dt / mass;

		x[i] = x[i] + xp[i] * dt;
		y[i] = y[i] + yp[i] * dt;

		/* z[i]  = F.real(); */
		/* zp[i] = F.imag(); */
		z[i] = F.real();
		zp[i] = -F.real()/GSL_CONST_MKSA_ELECTRON_CHARGE;
	}
	return 0;
}

int Ions::push_field(Field_Data &field, int z_step)
{
	double dt = _simparams->dt();
	/* Field_Interp fieldinterp(field, *gsl_interp2d_bicubic); */
	Field_Interp fieldinterp(field, *gsl_interp2d_bilinear);
	double Fx, Fy, Ex, Ey;
	int x_ind, y_ind;

	for (int i=0; i < n_pts-1; i++)
	{
		Ex = fieldinterp.Ex(x[i], y[i], z_step);
		Ey = fieldinterp.Ey(x[i], y[i], z_step);

		/* x_ind = round(x[i] / field.dxdi) + ( (field.x_pts-1)/2 ); */
		/* y_ind = round(y[i] / field.dydj) + ( (field.y_pts-1)/2 ); */
		/* Ex = field.Ex_ind(x_ind, y_ind, z_step); */
		/* Ey = field.Ey_ind(x_ind, y_ind, z_step); */

		Fx = -GSL_CONST_MKSA_ELECTRON_CHARGE * Ex;
		Fy = -GSL_CONST_MKSA_ELECTRON_CHARGE * Ey;

		xp[i] = xp[i] + Fx * dt / mass;
		yp[i] = yp[i] + Fy * dt / mass;

		x[i] = x[i] + xp[i] * dt;
		y[i] = y[i] + yp[i] * dt;

		/* z[i]  = x_ind; */
		/* zp[i] = y_ind; */

		z[i]  = Fx;
		zp[i] = Ex;
	}
	/* printf("Updated\n"); */
	return 0;
}
