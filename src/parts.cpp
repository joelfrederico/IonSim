#include "parts.h"
#include "simparams.h"
#include <exception>
#include "consts.h"
#include "support_func.h"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_histogram2d.h>
#include "scalar_data.h"

// ==============================
// Constructors
// ==============================
int Parts::_init()
{
	switch (type)
	{
		case PARTS_ION:
			_particle_charge = GSL_CONST_MKSA_FARADAY / GSL_CONST_NUM_AVOGADRO;
			break;
		case PARTS_E:
			_particle_charge = -GSL_CONST_MKSA_FARADAY / GSL_CONST_NUM_AVOGADRO;
			break;
	}

	return 0;
}

Parts::Parts(double _mass, long long _n_pts, parttype_t _type) : 
	mass(_mass),
	n_pts(_n_pts),
	type(_type)
{
	_init();
}

long long calc_n_pts(const SimParams &simparams, const parttype_t type)
{
	switch (type)
	{
		case PARTS_ION:
			return simparams.n_ions_node();
			break;
		case PARTS_E:
			return simparams.n_e_node();
			break;
	}
	return 0;
}

Parts::Parts(const SimParams &simparams, const parttype_t _type) : type(_type), mass(simparams.ion_mass()), n_pts(calc_n_pts(simparams, _type))
{
	_init();
}

ScalarData Parts::get_rho_dz(double z0, double z1, const SimParams &simparams)
{
	int x_pts = simparams.n_field_x;
	int y_pts = simparams.n_field_y;
	long double x_edge_mag = simparams.field_trans_wind;
	long double y_edge_mag = simparams.field_trans_wind;

	ScalarData rho(x_pts, y_pts, 1, x_edge_mag, y_edge_mag, 0);
	for (long long i=0; i < n_pts; i++)
	{
		if ((z0 <= z[i]) && (z[i] < z1))
		{
			rho.ind(rho.lt_x_ind(x[i]), rho.lt_y_ind(y[i]), rho.lt_z_ind(z[i]))++;
		}
	}

	rho *= _particle_charge;

	return rho;
}
