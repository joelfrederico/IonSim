#ifndef __SUPPORT_FUNC_H_INCLUDED__
#define __SUPPORT_FUNC_H_INCLUDED__

#include <mpi.h>
#include <string>
#include "parts.h"
#include <hdf5.h>
#include "field_data.h"
#include "field_comm.h"

namespace ionsim
{
	// ==================================
	// Methods
	// ==================================
	double GeV2gamma(double GeV);
	double gamma2GeV(double gamma);
	double gaussian(double mean, double sigma);

	double sr(double emit_n, double E, double n_p_cgs, double m_ion_amu);
	double nb_0(double q_tot, double sz, double sr);
	double nb_0(double q_tot, double sz, double emit_n, double E, double n_p_cgs, double m_ion_amu);
}

#endif
