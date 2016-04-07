#include "consts.h"
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

}
