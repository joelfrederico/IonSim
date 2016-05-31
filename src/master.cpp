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

int master(bool verbose)
{
	// ==============================
	// Starting gun
	// ==============================
	LoopComm loopcomm;
	MPI_Status status;

	if (verbose) printf("Master says: I am the MASTER!\n");
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
	Field_Data *field_interp;
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
	for (int step=0; step < simparams.n_steps; step++)
	{
		// ==============================
		// Allocate for this loop
		// ==============================
		field = new Field_Data(simparams);
		printf("Step: %d\n", step);

		// ==============================
		// Get fields from slaves
		// ==============================
		loopcomm.instruct(LOOP_GET_EFIELD);
		fieldcomm.recv_field_others_add(*field);

		// ==============================
		// Write electrons
		// ==============================
		loopcomm.instruct(LOOP_DUMP_E);
		loopcomm.send_slaves(step);

		// ==============================
		// Write total field
		// ==============================
		writer_s = new WriterSerial(simparams.filename);
		writer_s->writedata(step, *field);
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
		loopcomm.instruct(LOOP_DUMP_IONS);
		loopcomm.send_slaves(step);
		loopcomm.send_slaves(0);
		for (int z_step=0; z_step < simparams.n_field_z; z_step++)
		{
			std::cout << "Ion step: " << z_step << std::endl;
			loopcomm.instruct(LOOP_PUSH_IONS);
			/* MPI_Barrier(MPI_COMM_WORLD); */
			loopcomm.send_slaves(step);

			loopcomm.instruct(LOOP_DUMP_IONS);
			loopcomm.send_slaves(step);
			loopcomm.send_slaves(z_step+1);

		}

		// ==============================
		// Deallocate for this loop
		// ==============================
		delete field;
	}

	loopcomm.instruct(LOOP_KILL);


	return 0;
}
