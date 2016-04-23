#include "beam.h"
#include <math.h>
#include <gsl/gsl_const_mksa.h>

// ==============================
// Cov
// ==============================
int Cov::_index(int i, int j)
{
	return i + 2*j;
}

Cov::Cov()
{
	for (int i=0; i < 2; i++)
	{
		for (int j=0; j < 2; j++)
		{
			_cov[_index(i, j)] = 0;
		}
	}
}

Cov::Cov(double cov[2][2])
{
	for (int i=0; i < 2; i++)
	{
		for (int j=0; j < 2; j++)
		{
			_cov[_index(i, j)] = cov[i][j];
		}
	}
}

Cov::Cov(double cov[4])
{
	for (int i=0; i < 4; i++)
	{
		_cov[i] = cov[i];
	}
}

Cov::Cov(double xx, double xy, double yx, double yy)
{
	Cov out;
	out(0, 0) = xx;
	out(0, 1) = xy;
	out(1, 0) = yx;
	out(1, 1) = yy;
}

double & Cov::operator()(int i, int j)
{
	return _cov[_index(i, j)];
}

// ==============================
// Beam
// ==============================
Beam::Beam() {}

Beam::Beam(double beta, double alpha, Emit emit)
{
	_beta  = beta;
	_alpha = alpha;
	_emit = emit;
}

double Beam::alpha() const
{
	return _alpha;
}

double Beam::beta() const
{
	return _beta;
}

double Beam::sigma() const
{
	return sqrt(_beta*_emit.emit());
}

Cov Beam::cov() const
{
	Cov out;

	out(0, 0) = beta();
	out(0, 1) = out(1, 0) = -alpha();
	out(1, 1) = (1.0 + pow(alpha(), 2))/beta();
	double buf = _emit.emit();
	for (int i=0; i < 2; i++) 
	{
		for (int j=0; j < 2; j++)
		{
			out(i, j) *= buf;
		}
	}

	return out;

	/* output[0][0] = beta(); */
	/* output[0][1] = output[1][0] = -alpha(); */
	/* output[1][1] = (1.0+pow(alpha(), 2))/beta(); */
	/* double buf = _emit.emit(); */
	/* for (int i=0; i < 2; i++) */ 
	/* { */
	/* 	for (int j=0; j < 2; j++) */
	/* 	{ */
	/* 		output[i][j] *= _emit.emit(); */
	/* 	} */
	/* } */
	/* return; */
}

// ==============================
// Match
// ==============================
Match::Match(Plasma plasma, double E, Emit emit)
{
	_plasma = plasma;
	_E = E;
	_emit = emit;
}

double Match::beta()
{
	return 1.0 / sqrt( _plasma.k_ion(_E));
}

double Match::alpha()
{
	return 0;
}


// ==================================
// Plasma
// ==================================
Plasma::Plasma() {}

Plasma::Plasma(double n_p_cgs, double ion_mass_amu)
{
	_n_p          = n_p_cgs * 1e6;
	_ion_mass_amu = ion_mass_amu;
}

double Plasma::n_p()
{
	return _n_p;
}

double Plasma::m()
{
	return _ion_mass_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
}

double Plasma::w_p()
{
	return sqrt( n_p() * pow(GSL_CONST_MKSA_ELECTRON_CHARGE, 2.0) / (GSL_CONST_MKSA_MASS_ELECTRON * GSL_CONST_MKSA_VACUUM_PERMITTIVITY) );
}

double Plasma::k_ion(double E)
{
	return n_p() * pow(GSL_CONST_MKSA_ELECTRON_CHARGE, 2.0) / (2.0 * E * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT * GSL_CONST_MKSA_VACUUM_PERMITTIVITY);
}

