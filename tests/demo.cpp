#include "ebeam.h"
#include "fields.h"
#include "simparams.h"
#include "support_func.h"
#include <iostream>
#include <math.h>
#include <stdio.h>

int main(int argv, char **argc)
{
	long n_e             = 1e4;
	long n_ions          = 1e4;
	double q_tot         = 2e10;
	double radius        = 2.440175e-7*10;
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
	
	Emit emit;
	emit.set_emit_n(emit_n, E);
	Plasma plas(n_p_cgs, m_ion_amu);
	Match mat(plas, E, emit);
	Beam x_beam(mat.beta(), mat.alpha(), emit);
	Beam y_beam(mat.beta(), mat.alpha(), emit);

	unsigned long int s = 1;
	Ebeam ebeam(simparams, x_beam, y_beam, s);

	ionsim::overwrite_file_serial(filename);

	Field field(simparams);

	field.dump_serial(filename, 0);

	ebeam.field(field);

	field.dump_serial(filename, 1);

	return 0;
}
