#include "consts.h"
#include "ebeam.h"
#include "emit.h"
#include "fftw_classes.h"
#include "field_data.h"
#include "ions.h"
#include "loop_comm.h"
#include "mpi.h"
#include "scalar_data_comm.h"
#include "support_func.h"
#include "writer_parallel.h"
#include "writer_serial.h"
#include <fftw3-mpi.h>
#include <gflags/gflags.h>
#include <gsl/gsl_const.h>
#include <iomanip>
#include <math.h>
#include <sstream>
#include "slave_loops.h"

DECLARE_bool(verbose);

int slave()
{
	// ==================================
	// Initialize Variables
	// ==================================
	// FFTW variables
	ptrdiff_t N0, N1;

	// Communication variables
	LoopComm loopcomm;
	int buf;

	// Simulation metadata
	SimParams simparams_temp;

	// Loop variables
	unsigned int step_buf, substep_buf;
	bool loop_alive;

	// Old-style fields
	Field_Data *field;
	Field_Data *ion_field;
	
	// IO variables

	// Particle variables
	Emit emit;

	/* int buf, step_buf, substep_buf; */

	// ==================================
	// Receive starting gun
	// ==================================
	MPI_Recv(&buf, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if (FLAGS_verbose) printf("Slave %d says: **DUM DUM DUM DUM**\n", loopcomm.id);

	MPI_Barrier(MPI_COMM_WORLD);

	// ==================================
	// Initialize FFTW
	// ==================================
	fftwl_mpi_broadcast_wisdom(MPI_COMM_WORLD);

	// ==================================
	// Receive simparams
	// ==================================
	simparams_temp.bcast_receive();
	const SimParams simparams = simparams_temp;

	// ==================================
	// Set up fields
	// ==================================
	field = new Field_Data(simparams);
	ion_field = new Field_Data(simparams.n_field_x, simparams.n_field_y, 1, simparams.field_trans_wind, simparams.field_trans_wind, 0);
	ScalarData<ldouble> rho(simparams);
	ScalarData<ldouble> psi(simparams);

	// ==================================
	// Generate beam
	// ==================================
	// Beam metadata
	emit.set_emit_n(simparams.emit_n, simparams.E);
	Plasma plas(simparams.n_p_cgs, simparams.m_ion_amu);
	Match mat(plas, simparams.E, emit);

	// Beam transverse
	Beam x_beam(mat.beta(), mat.alpha(), emit);
	Beam y_beam(mat.beta(), mat.alpha(), emit);

	// Beam particles
	Ebeam ebeam(simparams, x_beam, y_beam);

	// ==================================
	// Generate ions
	// ==================================
	Ions ions(&simparams, plas);

	// ==================================
	// Slave loop
	// ==================================
	loop_alive = true;
	do
	{
		loopcomm.instruct(&buf);
		switch (buf)
		{
			// ==================================
			// Set electron loop iteration
			// ==================================
			case LOOP_START_E_ITER:
				loopcomm.recv_master(&step_buf);

				break;

			// ==================================
			// Set ion loop iteration
			// ==================================
			case LOOP_START_I_ITER:
				loopcomm.recv_master(&substep_buf);
				break;

			// ==================================
			// Terminate Loop
			// ==================================
			case LOOP_KILL:
				loop_alive = false;
				break;

			// ==================================
			// Write ions to file
			// ==================================
			case LOOP_DUMP_IONS:
				SL_dump_ions(simparams, loopcomm, step_buf, substep_buf, ions);

				break;

			// ==================================
			// Write electrons to file
			// ==================================
			case LOOP_DUMP_E:
				SL_dump_electrons(simparams, loopcomm, step_buf, ebeam);

				break;

			// ==================================
			// Push ions
			// ==================================
			case LOOP_PUSH_IONS:
				SL_push_ions(simparams, substep_buf, field, ions);
				break;

			// ==================================
			// Send field to Master
			// ==================================
			case LOOP_GET_EFIELD:
				SL_get_efield(simparams, ebeam, field);

				break;

			// ==================================
			// Retrieve field from Master
			// ==================================
			case LOOP_SEND_EFIELD:
				SL_send_efield(simparams, field);

				break;

			// ==================================
			// Retrieve current ion field
			// ==================================
			case LOOP_GET_IFIELD:
				SL_get_ifield(ion_field);

				break;

			case LOOP_RESET_IFIELD:
				delete ion_field;
				ion_field = new Field_Data(simparams);

				break;

			case LOOP_GET_RHO:
				SL_get_rho(step_buf, substep_buf, simparams, rho, ebeam);

				break;

			case LOOP_GET_FIELDS:
				// ==================================
				// Prepare for FFT
				// ==================================
				JTF_PRINT(Starting sub);
				N0 = N1 = 256;
				psifftw_base(loopcomm, N0, N1);

				JTF_PRINT_NOEND(Sub done) << " (id: " << loopcomm.id << ")" << std::endl;

				break;
		}

	} while ( loop_alive == true );

	// ==================================
	// Clean up
	// ==================================
	delete field;

	return 0;
}
