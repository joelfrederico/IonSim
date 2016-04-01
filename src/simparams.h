#ifndef __SIMPARAMS_H_INCLUDED__
#define __SIMPARAMS_H_INCLUDED__

#include <string>
#include <mpi.h>
#include "consts.h"

class SimParams
{
	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
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
			pushmethod_t _pushmethod,
			long _n_e,
			long _n_field_x,
			long _n_field_y,
			long _n_field_z,
			long _n_ions,
			std::string _filename
			);

		SimParams();

		// ==================================
		// Communications
		// ==================================
		int bcast_send() const;
		int bcast_receive();

		/* int write_attributes_parallel(Writer &writer) const; */

		template <class T>
		int bcast_send_wrap(T send) const
		{
			T buf = send;
			if (typeid(send) == typeid(double)) {
				MPI::COMM_WORLD.Bcast(&buf, 1, MPI::DOUBLE, 0);
			} else if (typeid(send) == typeid(int)) {
				MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
			} else if (typeid(send) == typeid(long)) {
				MPI::COMM_WORLD.Bcast(&buf, 1, MPI::LONG, 0);
			} else return 0;
			return 1;
		}

		// ==================================
		// Data
		// ==================================
		double E;
		double dt;
		double emit_n;
		double length;
		double m_ion_amu;
		double n_p_cgs;
		double q_tot;
		double radius;
		double sdelta;
		double sz;
		double t_tot;
		int n_steps;
		pushmethod_t pushmethod;
		long n_e;
		long n_field_x;
		long n_field_y;
		long n_field_z;
		long n_ions;
		std::string filename;
		double gamma_rel;

		// ==================================
		// Data methods
		// ==================================
		int z_cov(double (&out)[2][2]);
		double ion_mass() const;
};

#endif
