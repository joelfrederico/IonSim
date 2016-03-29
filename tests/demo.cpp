#include "ebeam.h"
#include "fields.h"
#include "simparams.h"
#include "support_func.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int main(int argv, char **argc)
{
	/* const int NMAX = 300; // maximum number of array points */
	double *x_array, *y_array, *z_array;
	double x_temp, y_temp, z_temp;

	double mag = 2;
	int npts = 201;
	double delta = double(npts)/(2.0*mag);
	double offset = (npts-1)/2;

	// Step 2a: declare the accelerator and the spline object.
	gsl_interp_accel *xaccel_ptr;
	gsl_interp_accel *yaccel_ptr;
	gsl_spline2d *spline_ptr;

	// Step 2b: allocate the accelerator and the spline object.
	xaccel_ptr = gsl_interp_accel_alloc ();
	yaccel_ptr = gsl_interp_accel_alloc ();

	spline_ptr = gsl_spline2d_alloc (gsl_interp2d_bicubic, npts, npts); // cubic spline

	// Interpolate y = sin(x^2) from 0 to 2 with 20 points

	x_array = new double[npts];
	y_array = new double[npts];
	z_array = new double[npts*npts];

	for (int i = 0; i < npts; i++)
	{
		x_temp = (i-offset) / delta;
		for (int j=0; j < npts; j++)
		{
			y_temp = (j-offset) / delta;

			x_array[i] = x_temp;
			y_array[j] = y_temp;

			z_temp = sin (x_temp * x_temp) + cos(y_temp*y_temp);

			gsl_spline2d_set(spline_ptr, z_array, i, j, z_temp);
		}
	}


	// Step 3: initialize the spline
	gsl_spline2d_init(spline_ptr, x_array, y_array, z_array, npts, npts);

	// Step 4: evaluate the spline and/or derivatives as needed
	double x        = 1; // test point
	double y        = -1; // test point
	int x_ind       = round(x*delta+offset);
	int y_ind       = round(y*delta+offset);
	double z        = gsl_spline2d_eval (spline_ptr, x, y, xaccel_ptr, yaccel_ptr);
	double z_deriv  = gsl_spline2d_eval_deriv_x (spline_ptr, x, y, xaccel_ptr, yaccel_ptr);
	double z_deriv2 = gsl_spline2d_eval_deriv_xx (spline_ptr, x, y, xaccel_ptr, yaccel_ptr);

	std::cout << "offset=" << offset << ", delta=" << delta << std::endl;
	std::cout << "x=" << x << ", y=" << y << ", z=" << z << ", z’=" << z_deriv << ", z’’=" << z_deriv2 << std::endl;
	std::cout << "x_set=" << x_ind << ", y_set=" << y_ind << ", z_set=" << z_array[y_ind*npts + x_ind] << std::endl;
	std::cout << "x_set=" << x_ind << ", y_set=" << y_ind << ", z_set=" << gsl_spline2d_get(spline_ptr, z_array, x_ind, y_ind) << std::endl;
	// Step 5: free the accelerator and spline object
	gsl_spline2d_free (spline_ptr);
	gsl_interp_accel_free (xaccel_ptr);
	gsl_interp_accel_free (yaccel_ptr);
	return (0); // successful completion
}

int oldmain(int argv, char **argc)
{
	double x_cov[2][2];
	double y_cov[2][2];

	long n_e             = 1e4;
	long n_ions          = 1e4;
	double q_tot         = 2e10;
	double radius        = 6.133e-9;
	double length        = 100e-6;
	double E             = 20.35;
	double emit_n        = 50e-6;
	double n_p_cgs       = 1e17;
	double m_ion_amu     = 1.00794;
	double sz            = 30e-6;
	double sdelta        = 0.04;
	double t_tot         = 1.58631e-12*2;
	int n_steps          = 100;
	double dt            = t_tot/n_steps;
	std::string filename = "output.h5";
	int runge_kutta      = 0;
	long n_field_x       = 101;
	long n_field_y       = 101;
	long n_field_z       = 101;
	unsigned long int s  = 1;

	int x_mag = 5;
	int y_mag = 1;

	Emit emit;
	emit.set_emit_n(emit_n, E);
	Plasma plas(n_p_cgs, m_ion_amu);
	Match mat(plas, E, emit);
	Beam x_beam(mat.beta(), mat.alpha(), emit);
	Beam y_beam(mat.beta(), mat.alpha(), emit);

	x_beam.cov(x_cov);
	y_beam.cov(y_cov);

	/* radius = x_beam.sigma()*10; */
	radius = 10;

	/* radius = x_cov[0][0]*0.00001; */
	printf("Radius: %e\n", radius);
	printf("X_cov: %e\n", x_cov[0][0]);

	SimParams simparams(
		E,
		dt,
		emit_n,
		length,
		m_ion_amu,
		n_p_cgs,
		q_tot,
		radius,
		sdelta,
		sz,
		t_tot,
		n_steps,
		runge_kutta,
		n_e,
		n_field_x,
		n_field_y,
		n_field_z,
		n_ions,
		filename
		);
	
	std::cout << simparams.radius << std::endl;

	/* Ebeam ebeam(simparams, x_cov[0][0]*x_mag*x_mag, y_cov[0][0]*y_mag*y_mag, s); */
	Ebeam ebeam(simparams, x_mag*x_mag, y_mag*y_mag, s);
	/* Ebeam ebeam(simparams, x_beam, y_beam, s); */

	ionsim::overwrite_file_serial(filename);

	Field field(simparams);

	std::cout << "Field x_edge_mag: " << field.x_edge_mag << std::endl;

	ebeam.field(field);
	field.dump_serial(filename, 0);

	return 0;
}
