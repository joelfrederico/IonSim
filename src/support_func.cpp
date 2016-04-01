#include "consts.h"
#include "field_data.h"
#include "support_func.h"
#include <hdf5.h>
#include <iomanip>
#include <math.h>
#include <mpi.h>
#include <sstream>
#include <string>

// Number of entries to write at once
// because HDF5 crashes if you try to
// do too much all at once.
const int MAX_N_WRITE = 1e5;

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

	int sendloop(const int &message)
	{
		int buf;
		buf = message;
		MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
		return 0;
	}

	int loop_get_fields(Field_Comm &fieldcomm, Field_Data &field)
	{
		sendloop(LOOP_GET_EFIELD);

		fieldcomm.recv_field_others(field);
		return 0;
	}

	int loop_push_ions(Field_Data &field)
	{
		sendloop(LOOP_PUSH_IONS);

		int p = MPI::COMM_WORLD.Get_size();
		for (int id=1; id < p; id++) {
			/* field.send_field(id); */
		}
		return 0;
	}

	int sendloop(const int &message, int step)
	{
		sendloop(message);
		MPI::COMM_WORLD.Bcast(&step, 1, MPI::INT, 0);
		return 0;
	}
}
