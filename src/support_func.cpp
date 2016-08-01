#include "consts.h"
#include "emit.h"
#include "field_data.h"
#include "support_func.h"
#include <hdf5.h>
#include <iomanip>
#include <math.h>
#include <mpi.h>
#include <sstream>
#include <string>


namespace ionsim
{
	// ==================================
	// Methods
	// ==================================
	double GeV2gamma(double GeV)
	{
		return GeV * 1e9 * GSL_CONST_MKSA_ELECTRON_VOLT / ELECTRON_REST_ENERGY;
	}

	double gamma2GeV(double gamma)
	{
		return gamma * ELECTRON_REST_ENERGY / (1e9 * GSL_CONST_MKSA_ELECTRON_VOLT);
	}

	double gaussian()
	{
		return 0;
	}

	double sr(double emit_n, double E, double n_p_cgs, double m_ion_amu)
	{
		Emit emit;
		emit.set_emit_n(emit_n, E);

		Plasma plas(n_p_cgs, m_ion_amu);
		Match mat(plas, E, emit);
		Beam x_beam(mat.beta(), mat.alpha(), emit);
		return x_beam.sigma();
	}

	double nb_0(double q_tot, double sz, double emit_n, double E, double n_p_cgs, double m_ion_amu)
	{
		double _sr = sr(emit_n, E, n_p_cgs, m_ion_amu);
		return nb_0(q_tot, sz, _sr);
	}

	double nb_0(double q_tot, double sz, double sr)
	{
		return q_tot / (pow(2*M_PI, 1.5) * sz * sr * sr);
	}
}
