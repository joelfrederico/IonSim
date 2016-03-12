#ifndef __BASECLASS_H_INCLUDED__
#define __BASECLASS_H_INCLUDED__

#include "consts.h"

class Parts
{
	protected:
		long _n_pts;
		double _mass;
		double_vec _x;
		double_vec _xp;
		double_vec _y;
		double_vec _yp;
		double_vec _z;
		double_vec _zp;
	public:
		Parts(long n_pts, double mass);

		long n_pts();
		const double_vec * x();
		const double_vec * xp();
		const double_vec * y();
		const double_vec * yp();
		const double_vec * z();
		const double_vec * zp();
};

#endif
