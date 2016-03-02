#include "baseclass.h"
#include <gsl/gsl_const_mksa.h>
#include "support_func.h"

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
Emit::Emit() 
{
}

void Emit::set_emit(double emit, double E_GeV)
{
	_emit = emit;
}

void Emit::set_emit_n(double emit_n, double E_GeV)
{
	_emit = emit_n / ionsim::GeV2gamma(E_GeV);
}

double Emit::emit()
{
	return _emit;
}
