#include "ebeam.h"
#include "fields.h"
#include "simparams.h"
#include "support_func.h"
#include <iostream>
#include <math.h>
#include <stdio.h>

int main(int argv, char **argc)
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
