#include "baseclass.h"
#include <gsl/gsl_const_mksa.h>
#include "support_func.h"
#include "consts.h"

// ==============================
// Parts
// ==============================
Parts::Parts(long n_pts, double mass)
{
	_n_pts = n_pts;
	_mass  = mass;

	_x.reserve(n_pts);
	_xp.reserve(n_pts);
	_y.reserve(n_pts);
	_yp.reserve(n_pts);
	_z.reserve(n_pts);
	_zp.reserve(n_pts);
}

long Parts::n_pts()
{
	return _n_pts;
}

const double_vec * Parts::x()
{
	return &_x;
}

const double_vec * Parts::xp()
{
	return &_xp;
}

const double_vec * Parts::y()
{
	return &_y;
}

const double_vec * Parts::yp()
{
	return &_yp;
}

const double_vec * Parts::z()
{
	return &_z;
}

const double_vec * Parts::zp()
{
	return &_zp;
}

