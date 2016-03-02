#include "baseclass.h"
#include <gsl/gsl_const_mksa.h>

// ==============================
// Parts
// ==============================
Parts::Parts(int n_pts)
{
	_n_pts = n_pts;
	_x.reserve(n_pts);
	_xp.reserve(n_pts);
	_y.reserve(n_pts);
	_yp.reserve(n_pts);
	_z.reserve(n_pts);
	_zp.reserve(n_pts);
}

// ==============================
// Emit
// ==============================
Emit::Emit() {}

Emit::Emit(double emit, double E, bool emit_n)
{
	if ( emit_n )
	{
		_emit = emit / (E * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT);
	} else {
		_emit = emit;
	}
}

double Emit::emit()
{
	return _emit;
}
