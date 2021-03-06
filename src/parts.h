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
		Parts(double _mass, unsigned long long _n_pts, parttype_t _type);
		Parts(const SimParams &simparams, parttype_t _type);

		// ==============================
		// Data members
		// ==============================
		const double mass;
		const unsigned long long n_pts;
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
		template<typename T>
		int get_rho_dz(const double z0, const double z1, ScalarData<T> &rho) const;

		template<typename T>
		std::vector<ScalarData<T>> get_J();
};

#endif
