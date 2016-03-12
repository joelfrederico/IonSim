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

class Plasma
{
	private:
		double _n_p;
		double _ion_mass_amu;
	public:
		Plasma();
		Plasma(double n_p_cgs, double ion_mass_amu);

		double n_p();
		double m();
		double w_p();
		double k_ion(double E);
};

class Match
{
	private:
		Plasma _plasma;
		double _E;
		Emit _emit;
	public:
		Match(Plasma plasma, double E, Emit emit);
		double beta();
		double alpha();
};

#endif
