#include <mpi.h>
#include "support_func.h"
#include "field_data.h"
#include "field_comm.h"
#include "writer_serial.h"
#include "loop_comm.h"
#include "consts.h"
#include "simparams.h"
#include <gflags/gflags.h>
#include "scalar_data_comm.h"
#include <fftw3-mpi.h>
#include "fftw_classes.h"
#include "mpi_vec.h"
#include "master_loop.h"

DECLARE_string(file);
DECLARE_bool(verbose);

int master()
{
	// ==============================
	// Initialize variables
	// ==============================
	// Communications variables
	LoopComm loopcomm;
	ScalarData_Comm scalarcomm;
	Field_Comm fieldcomm;

	// FFTW vars
	const char *wisdom_file = ".fftw-wisdom";

	// ==============================
	// Starting gun
	// ==============================
	if (FLAGS_verbose) printf("Master says: I am the MASTER!\n");
	for (int slave_id=1; slave_id < loopcomm.p; slave_id++)
	{
		MPI_Send(&slave_id, 1, MPI_INT, slave_id, slave_id*2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// ==============================
	// Load wisdom
	// ==============================
	if (fftwl_import_wisdom_from_filename(wisdom_file))
	{
		std::cout << "Wisdom imported successfully" << std::endl;
	} else {
		std::cout << "Wisdom not imported" << std::endl;
	}
	fftwl_mpi_broadcast_wisdom(MPI_COMM_WORLD);


	// ==============================
	// Load from file
	// ==============================
	SimParams *simparams_try;
	try
	{
		simparams_try = new SimParams(FLAGS_file);
	} catch (...) {
		loopcomm.instruct(LOOP_KILL);
		throw;
	}

	const SimParams simparams = *simparams_try;
	delete simparams_try;

	// ==============================
	// Initialize fields based on simparams
	// ==============================
	ScalarData<ldouble> rho(simparams);

	// ==============================
	// Send simparams everywhere
	// ==============================
	simparams.bcast_send();

	// ==============================
	// Overwrite current file
	// ==============================
	ML_overwrite_file(simparams);

	// ==============================
	// Loop over electron evolution
	// ==============================
	for (unsigned int e_step=0; e_step < simparams.n_steps; e_step++)
	{
		// ==============================
		// Allocate for this loop
		// ==============================
		Field_Data field(simparams);
		printf("Step: %d\n", e_step);

		// ==============================
		// Setup slave for e loop
		// ==============================
		loopcomm.instruct(LOOP_START_E_ITER);
		loopcomm.send_slaves(e_step);

		// ==============================
		// Get fields from slaves
		// ==============================
		loopcomm.instruct(LOOP_GET_EFIELD);
		fieldcomm.recv_field_others_add(field);

		// ==============================
		// Write electrons
		// ==============================
		loopcomm.instruct(LOOP_DUMP_E);

		// ==============================
		// Write total field
		// ==============================
		ML_write_field(simparams, e_step, field);

		// ==============================
		// Push field to slaves
		// ==============================
		loopcomm.instruct(LOOP_SEND_EFIELD);
		for (int id=1; id < loopcomm.p; id++)
		{
			fieldcomm.send_field(field, id);
		}

		// ==============================
		// Get rho
		// ==============================
		loopcomm.instruct(LOOP_GET_RHO);
		rho = 0;
		/* scalarcomm.recv_scalar_others_add(rho); */

		break;

		// ==============================
		// Get fields
		// ==============================

		/*
		// ==============================
		// Integrate ion motion
		// ==============================
		for (int z_step=0; z_step < simparams.n_field_z; z_step++)
		{
			// ==============================
			// Setup slave for e loop
			// ==============================
			loopcomm.instruct(LOOP_START_E_ITER);
			loopcomm.send_slaves(z_step);

			// ==============================
			// Record ions
			// ==============================
			loopcomm.instruct(LOOP_DUMP_IONS);

			// ==============================
			// Get current ion field
			// ==============================
			loopcomm.instruct(LOOP_GET_IFIELD);
			loopcomm.recv_ion_field_others_add(*ion_field);

			// Only push ions if flag is set and not on last step
			// (no point to push ions if not going to record)
			if ( (simparams.push_ions) && (z_step < (simparams.n_field_z-1)) )
			{
				std::cout << "Ion step: " << z_step << std::endl;
				loopcomm.instruct(LOOP_PUSH_IONS);
			}

		}
		*/
	}

	loopcomm.instruct(LOOP_KILL);


	return 0;
}
