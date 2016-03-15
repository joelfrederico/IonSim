#include "baseclass.h"
#include <exception>
#include "consts.h"
#include "support_func.h"
#include <gsl/gsl_const_mksa.h>

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
	n_ions(_n_ions),
	filename(_filename)
{}

int SimParams::z_cov(double (&out)[2][2])
{
	out[0][0] = pow(sz, 2);
	out[0][1] = out[1][0] = 0;
	out[1][1] = pow(sdelta, 2);
	return 0;
}

double SimParams::ion_mass()
{
	return m_ion_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
}

// ==============================
// Parts
// ==============================
Parts::Parts(SimParams simparams, const parttype _type) : type(_type), mass(simparams.ion_mass()), field(simparams.n_field_x, simparams.n_field_y)
{
	_simparams = &simparams;

	switch (type)
	{
		case ionsim::PARTS_ION:
			_n_pts = (*_simparams).n_e;
			break;
		case ionsim::PARTS_E:
			_n_pts = (*_simparams).n_ions;
			break;
	}

	x.reserve(_n_pts);
	xp.reserve(_n_pts);
	y.reserve(_n_pts);
	yp.reserve(_n_pts);
	z.reserve(_n_pts);
	zp.reserve(_n_pts);
}

long Parts::n_pts() const
{
	return _n_pts;
}

