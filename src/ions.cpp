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
#include "mpi.h"

// ===================================
// Ions
// ===================================
Ions::Ions(const SimParams *simparams, Plasma &plasma) : Parts(*simparams, PARTS_ION), _simparams(simparams)
{
	// ===================================
	// Init variables
	// ===================================
	bool keep_looking;
	const gsl_rng_type * T;
	gsl_rng * r;
	int s;
	int i_start          = 0;
	bool custom_particle = false;

	// ===================================
	// Store config info
	// ===================================
	_plasma      = &plasma;
	_radius      = simparams->radius;

	// ===================================
	// Set up random number generator env
	// ===================================
	// Reads env vars and sets lib vars
	gsl_rng_env_setup();

	// ===================================
	// Set up random number generator
	// ===================================
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	MPI_Comm_rank(MPI_COMM_WORLD, &s);
	gsl_rng_set(r, s);

	// ===================================
	// Inject custom particle for
	// debugging
	// ===================================
	i_start = 0;
	custom_particle = true;
	if (custom_particle)
	{
		x[0]  = 2e-6;
		xp[0] = 0;
		y[0]  = 0;
		yp[0] = 0;
		z[0]  = 0;
		zp[0] = 0;

		i_start = 1;
	}

	long long n_ions_node = _simparams->n_ions_node();

	// ===================================
	// Generate ions
	// ===================================
	for (int i=i_start; i < n_ions_node; i++)
	{
		/* x[i]     = gsl_ran_flat(r, -_radius, _radius); */
		/* y[i]     = gsl_ran_flat(r, -_radius, _radius); */

		// ===================================
		// Only accept ions inside radius
		// ===================================
		keep_looking = true;
		while (keep_looking)
		{
			x[i]     = gsl_ran_flat(r, -_radius, _radius);
			y[i]     = gsl_ran_flat(r, -_radius, _radius);
			if ((x[i]*x[i] + y[i]*y[i]) < (_radius*_radius))
			{
				keep_looking = false;
			}
		}

		// ===================================
		// All motion to zero; z meaningless
		// ===================================
		z[i] = 0;
		xp[i] = 0;
		yp[i] = 0;
		zp[i] = 0;
	}

	gsl_rng_free(r);
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
	const double dt = _simparams->dt();
	const double dz = _simparams->dz();
	/* Field_Interp fieldinterp(field, *gsl_interp2d_bicubic); */
	Field_Interp fieldinterp(field, *gsl_interp2d_bilinear);
	double Fx, Fy, Ex, Ey;
	int x_ind, y_ind;

	for (int i=0; i < n_pts-1; i++)
	{
		Ex = fieldinterp.Ex(x[i], y[i], z_step);
		Ey = fieldinterp.Ey(x[i], y[i], z_step);

		Fx = -GSL_CONST_MKSA_ELECTRON_CHARGE * Ex;
		Fy = -GSL_CONST_MKSA_ELECTRON_CHARGE * Ey;

		xp[i] = xp[i] + Fx * dt / mass;
		yp[i] = yp[i] + Fy * dt / mass;

		x[i] = x[i] + xp[i] * dt;
		y[i] = y[i] + yp[i] * dt;

		z[i]  = dz * (z_step+1);
		zp[i] = Fx;
	}

	return 0;
}


int Ions::field_Coulomb_sliced(Field_Data &field, int step)
{
	double dx, dy, dz;
	double drsq, dr;
	double x_e, y_e, z_e;
	double k_double;
	double sr_m;
	double srsq_macro;
	int k_start, k_end;

	sr_m = sr_macro();
	srsq_macro = sr_m*sr_m;

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	/* const double common   = _simparams.qpp_e()*GSL_CONST_MKSA_ELECTRON_CHARGE / (4*M_PI*dz*srsq_macro); */

	/* const double common_E = common / GSL_CONST_MKSA_VACUUM_PERMITTIVITY; */
	/* const double common_B = common * GSL_CONST_MKSA_SPEED_OF_LIGHT * GSL_CONST_MKSA_VACUUM_PERMEABILITY; */

	double temp_tran;
	double temp_tran_E;
	double temp_tran_B;

	for (int n=0; n < n_pts; n++) {
		x_e = x[n];
		y_e = y[n];
		z_e = z[n];

		if (_simparams->ion_z_bool)
		{
			k_start = 0;
			k_end = field.z_pts;
		} else {
			k_start = step;
			k_end = step+1;
		}

		for (int i=0; i < field.x_pts; i++)
		{
			dx = field.x_grid[i] - x_e;
			for (int j=0; j < field.y_pts; j++)
			{
				dy = field.y_grid[j] - y_e;
				for (int k=k_start; k < k_end; k++)
				{
					dz = field.z_grid[k] - z_e;
					drsq = dx*dx + dy*dy + dz*dz;
					dr   = sqrt(drsq);

					temp_tran = gsl_sf_exprel(- drsq / (2*srsq_macro));

					temp_tran_E = temp_tran * common_E;
					temp_tran_B = temp_tran * common_B;

					field.Ex_ind(i, j, k) += temp_tran_E * dx;
					field.Ey_ind(i, j, k) += temp_tran_E * dy;

					/* field.Bx_ind(i, j, k) += -temp_tran_B * dy; */
					/* field.By_ind(i, j, k) +=  temp_tran_B * dx; */
				}
			}
		}
	}


	return 0;
}

double np() const
{
	return _simparams->n_p_cgs * 
}

// ==================================
// Smoothing- see thesis
// ==================================
double Ions::sr_macro() const
{
	return 0.2347535410605456 / sqrt(n_0() * _simparams.dz());
}
