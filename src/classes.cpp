#include "classes.h"
#include "consts.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>

Beam::Beam(float beta, float alpha)
{
	_beta  = beta;
	_alpha = alpha;
}

float Beam::alpha()
{
	return _alpha;
}

float Beam::beta()
{
	return _beta;
}

Ebeam::Ebeam(int nparts, float q_tot, float E, float sig_delta, float beta, float alpha)
{
	lies = 1;
}

Ebeam::~Ebeam()
{
	printf("%s", "lies");
}

Emit::Emit() {}

Emit::Emit(float emit, float E, bool emit_n)
{
	if ( emit_n )
	{
		_emit = emit / (E * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT);
	} else {
		_emit = emit;
	}
}

Plasma::Plasma() {}

Plasma::Plasma(float n_p_cgs, float ion_mass_amu)
{
	_n_p = n_p_cgs * 1e6;
	_ion_mass_amu = ion_mass_amu;
}

float Plasma::n_p()
{
	return _n_p;
}

float Plasma::m()
{
	return _ion_mass_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
}

float Plasma::w_p()
{
	return sqrt( n_p() * pow(GSL_CONST_MKSA_ELECTRON_CHARGE, 2.0) / (GSL_CONST_MKSA_MASS_ELECTRON * GSL_CONST_MKSA_VACUUM_PERMITTIVITY) );
}

float Plasma::k_ion(float E)
{
	return n_p() * pow(GSL_CONST_MKSA_ELECTRON_CHARGE, 2.0) / (2.0 * E * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT * GSL_CONST_MKSA_VACUUM_PERMITTIVITY);
}

Match::Match(Plasma plasma, float E, Emit emit)
{
	_plasma = plasma;
	_E = E;
	_emit = emit;
}

float Match::beta()
{
	return 1.0 / sqrt( _plasma.k_ion(_E));
}
