#ifndef __IONS_H_INCLUDED__
#define __IONS_H_INCLUDED__

#include "parts.h"

#include "mpi.h"
#include "beam.h"
#include "field_data.h"

class Ions : public Parts
{
	private:
		double _radius;
		double _part_charge;
		const SimParams *_simparams;

		Plasma * _plasma;

	public:
		Ions(const SimParams *simparams, Plasma &plasma, int n_pts, double radius, double length);

		int push(double nb_0, double sig_r);
		int push_simple(double nb_0, double sig_r);
		int push_field(Field_Data &field, int z_step);
};

#endif
