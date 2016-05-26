#include "support_func.h"
#include "support.h"

SimParams simparams_gen()
{
	std::string filename;
	// ==============================
	// Generate beam
	// ==============================
	long n_e                = 1e6;
	long n_ions             = 1e4;
	double q_tot            = 2e10;
	double radius           = 2.4276628847185805e-06;
	double sz               = 30e-6;
	double length           = 100e-6;
	double E                = 20.35;
	double emit_n           = 50e-6;
	double n_p_cgs          = 1e17;
	double m_ion_amu        = 1.00794;
	double sdelta           = 0.04;
	zdist_t zdist           = Z_DIST_FLAT;
	int n_steps             = 1;
	pushmethod_t pushmethod = PUSH_SIMPLE;

	if (pushmethod == PUSH_SIMPLE)
	{
		filename = "simple.h5";
	} else if (pushmethod == PUSH_FIELD) {
		filename = "field.h5";
	} else {
		filename = "output.h5";
	}

	long n_field_x          = 101;
	long n_field_y          = 101;
	long n_field_z          = 51;
	double field_trans_wind = radius;

	double sr = ionsim::sr(emit_n, E, n_p_cgs, m_ion_amu);
	double nb_0 = ionsim::nb_0(q_tot, sz, sr);

	double z_end = (11.1367*GSL_CONST_MKSA_SPEED_OF_LIGHT / GSL_CONST_MKSA_ELECTRON_CHARGE) * sqrt(GSL_CONST_MKSA_VACUUM_PERMITTIVITY * m_ion_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS / nb_0);

	q_tot *= z_end / sz;
	sz = z_end;

	return SimParams(
		E,
		emit_n,
		length,
		m_ion_amu,
		n_p_cgs,
		q_tot,
		radius,
		sz,
		sdelta,
		zdist,
		n_steps,
		pushmethod,
		n_e,
		n_field_x,
		n_field_y,
		n_field_z,
		field_trans_wind,
		z_end,
		n_ions,
		filename
		);
}
