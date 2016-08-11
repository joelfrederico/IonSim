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

DECLARE_string(file);
DECLARE_bool(verbose);

int master()
{
	// ==============================
	// Starting gun
	// ==============================
	const char *wisdom_file = ".fftw-wisdom";
	LoopComm loopcomm;
	ScalarData_Comm scalarcomm;

	if (FLAGS_verbose) printf("Master says: I am the MASTER!\n");
	for (int slave_id=1; slave_id < loopcomm.p; slave_id++)
	{
		MPI_Send(&slave_id, 1, MPI_INT, slave_id, slave_id*2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	// ==============================
	// Set up sim
	// ==============================
	int max_sim_steps = 100;
	int max_ion_steps = 100;
	
	// ==============================
	// Load from file
	// ==============================
	const SimParams simparams(FLAGS_file);

	// ==============================
	// Initialize fields
	// ==============================
	long long rho_size;
	ScalarData rho(simparams);
	ScalarData psi(simparams);
	Field_Data *field;
	Field_Comm fieldcomm;
	long long local_n0, local_0_start;
	long double *buf;

	// ==============================
	// Try to read FFT wisdom
	// ==============================
	fftwl_mpi_init();
	if (fftwl_import_wisdom_from_filename(".fftw_wisdom"))
	{
		std::cout << "Wisdom imported successfully" << std::endl;
	} else {
		std::cout << "Wisdom not imported" << std::endl;
	}

	fftwl_mpi_broadcast_wisdom(MPI_COMM_WORLD);

	// ==============================
	// Send simparams everywhere
	// ==============================
	simparams.bcast_send();

	// ==============================
	// Overwrite current file
	// ==============================
	WriterSerial *writer_s;
	writer_s = new WriterSerial(simparams.filename, true);
	(*writer_s).write_attributes(simparams);
	delete writer_s;

	// ==============================
	// Loop over electron evolution
	// ==============================
	for (int e_step=0; e_step < simparams.n_steps; e_step++)
	{
		// ==============================
		// Allocate for this loop
		// ==============================
		field     = new Field_Data(simparams);
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
		fieldcomm.recv_field_others_add(*field);

		// ==============================
		// Write electrons
		// ==============================
		loopcomm.instruct(LOOP_DUMP_E);

		// ==============================
		// Write total field
		// ==============================
		writer_s = new WriterSerial(simparams.filename);
		writer_s->writedata(e_step, *field);
		delete writer_s;

		// ==============================
		// Push field to slaves
		// ==============================
		loopcomm.instruct(LOOP_SEND_EFIELD);
		for (int id=1; id < loopcomm.p; id++)
		{
			fieldcomm.send_field(*field, id);
		}

		// ==============================
		// Get rho
		// ==============================
		loopcomm.instruct(LOOP_GET_RHO);
		rho = 0;
		scalarcomm.recv_scalar_others_add(rho);

		// ==============================
		// Get fields
		// ==============================
		loopcomm.instruct(LOOP_GET_FIELDS);
		for (int id=1; id < loopcomm.p; id++)
		{
			// Load wisdom and distribute
			if (fftwl_import_wisdom_from_filename(wisdom_file))
			{
				std::cout << "Wisdom import succeeded" << std::endl;
			} else {
				std::cout << "Wisdom import failed" << std::endl;
			}

			fftwl_recv_local_size(local_n0, local_0_start, id);

			rho_size = local_n0*rho.y_pts;

			// Get pointer to proper location in array
			buf = (rho.data.data() + local_0_start);

			// Send rho
			MPI_Send(&buf, local_n0*rho.y_pts, MPI_LONG_DOUBLE, id, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);
		}
		
		// Let slaves do FFT, then collect psi
		for (int id=1; id < loopcomm.p; id++)
		{
			fftwl_recv_local_size(local_n0, local_0_start, id);
			rho_size = local_n0*rho.y_pts;

			// Receive psi
			buf = (psi.data.data() + local_0_start);
			MPI_Recv(&buf, rho_size, MPI_LONG_DOUBLE, id, TAG_LOOP_MESSAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		}
		// Receive wisdom and save
		fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);
		fftwl_export_wisdom_to_filename(wisdom_file);

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

		// ==============================
		// Deallocate for this loop
		// ==============================
		delete field;
	}

	loopcomm.instruct(LOOP_KILL);


	return 0;
}
