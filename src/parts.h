#ifndef __PARTS_H_INCLUDED__
#define __PARTS_H_INCLUDED__

#include "simparams.h"
#include "consts.h"

class Parts
{
	private:
		int _vec_reserve();
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

		double_vec x;
		double_vec xp;
		double_vec y;
		double_vec yp;
		double_vec z;
		double_vec zp;

};

#endif
