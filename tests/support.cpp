#include "support_func.h"
#include "support.h"

SimParams simparams_gen()
{
	SimParams simparams;
	std::string filename;
	// ==============================
	// Generate beam
	// ==============================
	simparams.n_e        = 1e6;
	simparams.n_ions     = 1e4;
	simparams.q_tot      = 2e10;
	simparams.radius     = 2.4276628847185805e-06;
	simparams.sz         = 30e-6;
	simparams.length     = 100e-6;
	simparams.E          = 20.35;
	simparams.emit_n     = 50e-6;
	simparams.n_p_cgs    = 1e17;
	simparams.m_ion_amu  = 1.00794;
	simparams.sdelta     = 0.04;
	simparams.zdist      = Z_DIST_FLAT;
	simparams.n_steps    = 1;
	simparams.pushmethod = PUSH_SIMPLE;

	if (simparams.pushmethod == PUSH_SIMPLE)
	{
		simparams.filename = "simple.h5";
	} else if (simparams.pushmethod == PUSH_FIELD) {
		simparams.filename = "field.h5";
	} else {
		simparams.filename = "output.h5";
	}

	simparams.n_field_x        = 101;
	simparams.n_field_y        = 101;
	simparams.n_field_z        = 51;
	simparams.field_trans_wind = simparams.radius;

	double sr   = ionsim::sr(simparams.emit_n, simparams.E, simparams.n_p_cgs, simparams.m_ion_amu);
	double nb_0 = ionsim::nb_0(simparams.q_tot, simparams.sz, sr);

	simparams.z_end = (11.1367*GSL_CONST_MKSA_SPEED_OF_LIGHT / GSL_CONST_MKSA_ELECTRON_CHARGE) * sqrt(GSL_CONST_MKSA_VACUUM_PERMITTIVITY * simparams.m_ion_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS / nb_0);

	simparams.q_tot *= simparams.z_end / simparams.sz;
	simparams.sz = simparams.z_end;

	return simparams;
}
