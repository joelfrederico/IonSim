#ifndef __BASECLASS_H_INCLUDED__
#define __BASECLASS_H_INCLUDED__

#include "consts.h"

class Emit
{
	private:
		double _emit;
	public:
		Emit();
		Emit(double emit, double E, bool emit_n);

		double emit();
};

class Parts
{
	protected:
		int _n_pts;
		double_vec _x;
		double_vec _xp;
		double_vec _y;
		double_vec _yp;
		double_vec _z;
		double_vec _zp;
	public:
		Parts(int n_pts);
};

#endif
