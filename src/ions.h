#ifndef __IONS_H_INCLUDED__
#define __IONS_H_INCLUDED__

#include "parts.h"

#include "mpi.h"
#include "beam.h"

class Ions : public Parts
{
	private:
		double _radius;
		double _part_charge;

		Plasma * _plasma;

	public:
		Ions(SimParams &simparams, Plasma &plasma, int n_pts, double radius, double length);

		int dump(std::string const &filename, int step, MPI::Intracomm &comm);
		int push(double dt, double nb_0, double sig_r);
		int push_simple(double dt, double nb_0, double sig_r);
};

#endif
