#ifndef __PARTS_H_INCLUDED__
#define __PARTS_H_INCLUDED__

#include "simparams.h"
#include "consts.h"
#include "scalar_data.h"
#include <vector>

class Parts
{
	private:
		double _particle_charge;

		int _init();

	public:
		// ==============================
		// Constructors
		// ==============================
		Parts(double _mass, long long _n_pts, parttype_t _type);
		Parts(const SimParams &simparams, parttype_t _type);

		// ==============================
		// Data members
		// ==============================
		const double mass;
		const long long n_pts;
		const parttype_t type;

		ldouble_vec x;
		ldouble_vec xp;
		ldouble_vec y;
		ldouble_vec yp;
		ldouble_vec z;
		ldouble_vec zp;

		// ==============================
		// Data methods
		// ==============================
		ScalarData get_rho_dz(double z0, double z1, const SimParams &simparams);
		std::vector<ScalarData> get_J();
};

#endif
