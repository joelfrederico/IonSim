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
	x.resize(n_pts);
	xp.resize(n_pts);
	y.resize(n_pts);
	yp.resize(n_pts);
	z.resize(n_pts);
	zp.resize(n_pts);

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

Parts::Parts(double _mass, unsigned long long _n_pts, parttype_t _type) : 
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

template<typename T>
int Parts::get_rho_dz(const double z0, const double z1, ScalarData<T> &rho, const SimParams &simparams) const
{
	int x_pts = simparams.n_field_x;
	int y_pts = simparams.n_field_y;
	long double x_edge_mag = simparams.field_trans_wind;
	long double y_edge_mag = simparams.field_trans_wind;
	int x_ind, y_ind, z_ind, e_x, e_y, e_z;

	e_z = rho.lt_z_ind_e((z1+z0)/2, z_ind);
	for (long long i=0; i < n_pts; i++)
	{
		if ((z0 <= z[i]) && (z[i] < z1))
		{
			e_x = rho.lt_x_ind_e(x[i], x_ind);
			e_y = rho.lt_y_ind_e(y[i], y_ind);
			if ((e_x==0) && (e_y==0))
			{
				rho.ind(x_ind, y_ind, z_ind)++;
			}
		}
	}

	rho *= _particle_charge;

	return 0;
}

template int Parts::get_rho_dz(const double z0, const double z1, ScalarData<ldouble> &rho, const SimParams &simparams) const;
