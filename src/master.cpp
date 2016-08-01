#include "mpi.h"
#include "support_func.h"
#include "field_data.h"
#include "field_comm.h"
#include "writer_serial.h"
#include "loop_comm.h"
#include "consts.h"
#include "simparams.h"
#include <gflags/gflags.h>

DECLARE_string(file);
DECLARE_bool(verbose);

int master()
{
	// ==============================
	// Starting gun
	// ==============================
	LoopComm loopcomm;
	MPI_Status status;

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
	Field_Data *field;
	Field_Data *ion_field;
	Field_Comm fieldcomm;

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
		ion_field = new Field_Data(simparams);
		printf("Step: %d\n", e_step);

		// ==============================
		// Get fields from slaves
		// ==============================
		loopcomm.instruct(LOOP_GET_EFIELD);
		fieldcomm.recv_field_others_add(*field);

		// ==============================
		// Write electrons
		// ==============================
		loopcomm.instruct(LOOP_DUMP_E);
		loopcomm.send_slaves(e_step);

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
		// Integrate ion motion
		// ==============================
		for (int z_step=0; z_step < simparams.n_field_z; z_step++)
		{
			// ==============================
			// Record ions
			// ==============================
			loopcomm.instruct(LOOP_DUMP_IONS);
			loopcomm.send_slaves(e_step);
			loopcomm.send_slaves(z_step);

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
				loopcomm.send_slaves(z_step);
			}

		}

		// ==============================
		// Deallocate for this loop
		// ==============================
		delete ion_field;
		delete field;
	}

	loopcomm.instruct(LOOP_KILL);


	return 0;
}
