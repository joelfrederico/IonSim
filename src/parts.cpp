#include "parts.h"
#include "simparams.h"
#include <exception>
#include "consts.h"
#include "support_func.h"
#include <gsl/gsl_const_mksa.h>

// ==============================
// Constructors
// ==============================
Parts::Parts(double _mass, long _n_pts, parttype _type) : mass(_mass), n_pts(_n_pts), type(_type)
{
	_vec_reserve();
}

long calc_n_pts(const SimParams &simparams, const parttype type)
{
	switch (type)
	{
		case ionsim::PARTS_ION:
			return simparams.n_e;
			break;
		case ionsim::PARTS_E:
			return simparams.n_ions;
			break;
	}
	return 0;
}

Parts::Parts(const SimParams &simparams, const parttype _type) : type(_type), mass(simparams.ion_mass()), n_pts(calc_n_pts(simparams, _type))
{
	_vec_reserve();
}

// ==============================
// Private methods
// ==============================
int Parts::_vec_reserve()
{
	// ==============================
	// Free unneeded space
	// ==============================
	if (x.size() > n_pts)
	{
		x.resize(n_pts);
		xp.resize(n_pts);
		y.resize(n_pts);
		yp.resize(n_pts);
		z.resize(n_pts);
		zp.resize(n_pts);
	}

	// ==============================
	// Reserve needed space
	// ==============================
	x.reserve(n_pts);
	xp.reserve(n_pts);
	y.reserve(n_pts);
	yp.reserve(n_pts);
	z.reserve(n_pts);
	zp.reserve(n_pts);

	return 0;
}

