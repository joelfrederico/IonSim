#ifndef __PARTS_H_INCLUDED__
#define __PARTS_H_INCLUDED__

#include "simparams.h"
#include "consts.h"

class Parts
{
	protected:
		SimParams *_simparams;
		long _n_pts;

	public:
		Parts(SimParams simparams, parttype _type);

		const parttype type;
		const double mass;

		double_vec x;
		double_vec xp;
		double_vec y;
		double_vec yp;
		double_vec z;
		double_vec zp;

		long n_pts() const;
};

#endif
