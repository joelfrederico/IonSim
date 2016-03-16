#include "baseclass.h"
#include "simparams.h"
#include <exception>
#include "consts.h"
#include "support_func.h"
#include <gsl/gsl_const_mksa.h>

// ==============================
// Parts
// ==============================
Parts::Parts(SimParams simparams, const parttype _type) : type(_type), mass(simparams.ion_mass())
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

