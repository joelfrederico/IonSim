#ifndef __BASECLASS_H_INCLUDED__
#define __BASECLASS_H_INCLUDED__

#include "consts.h"
#include "fields.h"

class SimParams
{
	protected:

	public:
		SimParams(
			double _E,
			double _dt,
			double _emit_n,
			double _length,
			double _m_ion_amu,
			double _n_p_cgs,
			double _q_tot,
			double _radius,
			double _sdelta,
			double _sz,
			double _t_tot,
			int _n_steps,
			int _runge_kutta,
			long _n_e,
			long _n_field_x,
			long _n_field_y,
			long _n_ions,
			std::string _filename
			);
		// SimParams(const SimParams &rhs);
		// SimParams & operator=(const SimParams &rhs);

		const double E;
		const double dt;
		const double emit_n;
		const double length;
		const double m_ion_amu;
		const double n_p_cgs;
		const double q_tot;
		const double radius;
		const double sdelta;
		const double sz;
		const double t_tot;
		const int n_steps;
		const int runge_kutta;
		const long n_e;
		const long n_field_x;
		const long n_field_y;
		const long n_ions;
		const std::string filename;

		int z_cov(double (&out)[2][2]);
		double ion_mass();
};

class Parts
{
	protected:
		SimParams *_simparams;
		long _n_pts;
		Field field;

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
