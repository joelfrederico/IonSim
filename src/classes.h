#ifndef __CLASSES_H_INCLUDED__
#define __CLASSES_H_INCLUDED__

#include <gsl/gsl_rng.h>

class Beam
{
	private:
		float _alpha;
		float _beta;
	public:
		Beam();
		Beam(float beta, float alpha);

		float alpha();
		float beta();
};

class Ebeam
{
	private:
		int lies;
	public:
		/// Create Ebeam class.
		Ebeam(int nparts, float q_tot, float E, float sig_delta, float beta, float alpha);

		~Ebeam();
};

class Emit
{
	private:
		float _emit;
	public:
		Emit();
		Emit(float emit, float E, bool emit_n);
};

class Plasma
{
	private:
		float _n_p;
		float _ion_mass_amu;
	public:
		Plasma();
		Plasma(float n_p_cgs, float ion_mass_amu);

		float n_p();
		float m();
		float w_p();
		float k_ion(float E);
};

class Match
{
	private:
		Plasma _plasma;
		float _E;
		Emit _emit;
	public:
		Match(Plasma plasma, float E, Emit emit);
		float beta();
};

#endif
