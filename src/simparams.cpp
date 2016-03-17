#include "simparams.h"
#include <gsl/gsl_const_mksa.h>
#include <math.h>

// ==================================
// Constructors, Destructor
// ==================================
SimParams::SimParams(
	double _E,
	double _dt,
	double _emit_n,
	double _length,
	double _m_ion_amu,
	double _n_p_cgs,
	double _q_tot,
	double _radius,
	double _sdelta,
	double _sz,
	double _t_tot,
	int _n_steps,
	int _runge_kutta,
	long _n_e,
	long _n_field_x,
	long _n_field_y,
	long _n_field_z,
	long _n_ions,
	std::string _filename
	) : 
	E(_E),
	dt(_dt),
	emit_n(_emit_n),
	length(_length),
	m_ion_amu(_m_ion_amu),
	n_p_cgs(_n_p_cgs),
	q_tot(_q_tot),
	radius(_radius),
	sdelta(_sdelta),
	sz(_sz),
	t_tot(_t_tot),
	n_steps(_n_steps),
	runge_kutta(_runge_kutta),
	n_e(_n_e),
	n_field_x(_n_field_x),
	n_field_y(_n_field_y),
	n_field_z(_n_field_z),
	n_ions(_n_ions),
	filename(_filename)
{}

// ==================================
// Public methods
// ==================================
int SimParams::z_cov(double (&out)[2][2])
{
	out[0][0] = pow(sz, 2);
	out[0][1] = out[1][0] = 0;
	out[1][1] = pow(sdelta, 2);
	return 0;
}

double SimParams::ion_mass() const
{
	return m_ion_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
}
