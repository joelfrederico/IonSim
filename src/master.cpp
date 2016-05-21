#include "mpi.h"
#include "support_func.h"
#include "field_data.h"
#include "field_comm.h"
#include "writer_serial.h"
#include "loop_comm.h"
#include "consts.h"

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
	// Generate beam
	// ==============================
	long n_e                = 1e4;
	long n_ions             = 1e5;
	double q_tot            = 2e10;
	double radius           = 2.4276628847185805e-06;
	double sz               = 30e-6;
	double length           = 100e-6;
	double E                = 20.35;
	double emit_n           = 50e-6;
	double n_p_cgs          = 1e17;
	double m_ion_amu        = 1.00794;
	double sdelta           = 0.04;
	zdist_t zdist           = Z_DIST_FLAT;
	int n_steps             = 1;
	std::string filename    = "output.h5";
	pushmethod_t pushmethod = PUSH_SIMPLE;
	/* pushmethod_t pushmethod = PUSH_FIELD; */
	long n_field_x          = 51;
	long n_field_y          = 51;
	long n_field_z          = 21;
	double field_trans_wind = radius * 10;

	double sr = ionsim::nb_0(q_tot, sz, emit_n, E, n_p_cgs, m_ion_amu);
	double nb_0 = ionsim::nb_0(q_tot, sz, sr);

	double z_end = (11.1367*GSL_CONST_MKSA_SPEED_OF_LIGHT / GSL_CONST_MKSA_ELECTRON_CHARGE) * sqrt(GSL_CONST_MKSA_VACUUM_PERMITTIVITY * m_ion_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS / nb_0);

	q_tot *= z_end / sz;
	sz = z_end;

	const SimParams simparams(
		E,
		emit_n,
		length,
		m_ion_amu,
		n_p_cgs,
		q_tot,
		radius,
		sz,
		sdelta,
		zdist,
		n_steps,
		pushmethod,
		n_e,
		n_field_x,
		n_field_y,
		n_field_z,
		field_trans_wind,
		z_end,
		n_ions,
		filename
		);

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
	writer_s = new WriterSerial(filename, true);
	(*writer_s).write_attributes(simparams);
	delete writer_s;

	// ==============================
	// Loop over electron evolution
	// ==============================
	for (int step=0; step < n_steps; step++)
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
		writer_s = new WriterSerial(filename);
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
		for (int z_step=0; z_step < n_field_z; z_step++)
		{
			loopcomm.instruct(LOOP_PUSH_IONS);
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
