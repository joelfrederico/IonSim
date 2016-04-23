#ifndef __BEAM_H_INCLUDED__
#define __BEAM_H_INCLUDED__

#include "emit.h"

class Cov
{
	private:
		double _cov[4];

		int _index(int i, int j);

	public:
		Cov();
		Cov(double cov[2][2]);
		Cov(double arr[4]);
		Cov(double xx, double xy, double yx, double yy);

		double & operator()(int i, int j);
};

class Beam
{
	private:
		double _alpha;
		double _beta;
		Emit _emit;
	public:
		Beam();
		Beam(double beta, double alpha, Emit emit);

		double alpha() const;
		double beta() const;
		double sigma() const;
		Cov cov() const;

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
