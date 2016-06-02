#ifndef __IONS_H_INCLUDED__
#define __IONS_H_INCLUDED__

#include "parts.h"

#include "mpi.h"
#include "beam.h"
#include "field_data.h"

std::complex<double> F_r(double x, double y, const SimParams &simparams);

class Ions : public Parts
{
	private:
		double _radius;
		const SimParams *_simparams;

		Plasma * _plasma;

	public:
		Ions(const SimParams *simparams, Plasma &plasma);

		int push(double nb_0, double sig_r);
		int push_simple(double nb_0, double sig_r);
		int push_field(Field_Data &field, int z_step);
};

#endif
