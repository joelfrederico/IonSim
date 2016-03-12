#ifndef __BEAM_H_INCLUDED__
#define __BEAM_H_INCLUDED__

#include "emit.h"

class Beam
{
	private:
		double _alpha;
		double _beta;
		Emit _emit;
	public:
		Beam();
		Beam(double beta, double alpha, Emit emit);

		double alpha();
		double beta();
		double sigma();
		void cov(double output[2][2]);

};

#endif
