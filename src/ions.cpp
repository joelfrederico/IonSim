#include "ions.h"
#include "fields.h"
#include <complex>
#include <sstream>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include "support_func.h"
#include "simparams.h"

// ==============================
// Ions
// ==============================
Ions::Ions(const SimParams &simparams, Plasma &plasma, int n_pts, double radius, double length) : Parts(simparams, ionsim::PARTS_ION)
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
	/* gsl_rng_set(r, MPI::COMM_WORLD.Get_rank() + 1); */

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

int Ions::dump_parallel(std::string const &filename, int step, MPI::Intracomm &comm)
{
	ionsim::dump_parallel(filename, step, "ions", comm, *this);
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
	double params[3] = {nb_0, sig_r, mass};
	double t = 0;
	int p = MPI::COMM_WORLD.Get_rank();

	gsl_odeiv2_system sys = {func, NULL, 4, &params};
	gsl_odeiv2_evolve *e  = gsl_odeiv2_evolve_alloc(4);
	gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, 4);

	/* gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, dt, sig_r/100, 1e-8); */
	
	for (int i=0; i < n_pts; i++)
	{
		double y[] = {x[i], y[i], xp[i], yp[i]};
		double yerr[4];
		/* printf("%d Trying %d\n", p, i); */
		/* gsl_odeiv2_driver_apply_fixed_step(d, &t, dt, 4, y); */
		gsl_odeiv2_step_apply(step, t, dt, y, yerr, NULL, NULL, &sys);
		x[i]  = y[0];
		xp[i] = y[2];
		y[i]  = y[1];
		yp[i] = y[3];

		F = F_r(x[i], y[i], nb_0, sig_r);
		z[i]  = std::abs(F);
		zp[i] = mass;
		/* printf("%d Got here\n", p); */
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
	}
	return 0;
}

int Ions::push_field(double dt, Field &field)
{
	// ==============================
	// Gather 
	// ==============================
	double Fx, Fy;
	/* std::complex<double> F; */
	for (int i=0; i < n_pts; i++)
	{
		Fx = GSL_CONST_MKSA_ELECTRON_CHARGE * field.Ex(x[i], y[i]);
		Fy = GSL_CONST_MKSA_ELECTRON_CHARGE * field.Ey(x[i], y[i]);

		x[i] += xp[i] * dt;
		y[i] += yp[i] * dt;

		xp[i] += Fx * dt / mass;
		yp[i] += Fy * dt / mass;
		/* if (i==0) { */
		/* 	printf("Fx = %e, i: %d, x: %e, y: %e\n", Fx, i, x[i], y[i]); */
		/* 	printf("Fx = %e\n", field.Ex_ind(0, 0)); */
		/* 	printf("x = %e\n", field.i_to_x(0)); */
		/* } */
	}
	/* printf("Updated\n"); */
	return 0;
}

