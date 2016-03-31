#include "mpi.h"
#include "support_func.h"
#include "field_data.h"
#include "field_comm.h"

int master(int &p, bool verbose)
{
	// ==============================
	// Starting gun
	// ==============================
	MPI::Status status;

	if (verbose) printf("Master says: I am the MASTER!\n");
	for (int slave_id=1; slave_id < p; slave_id++)
	{
		MPI::COMM_WORLD.Send(&slave_id, 1, MPI::INT, slave_id, slave_id*2);
	}
	// ==============================
	// Set up sim
	// ==============================
	int max_sim_steps = 100;
	int max_ion_steps = 100;
	
	// ==============================
	// Generate beam
	// ==============================
	
	long n_e                = 1e5;
	long n_ions             = 1e4;
	double q_tot            = 2e10;
	double radius           = 2.4276628847185805e-06;
	double length           = 100e-6;
	double E                = 20.35;
	double emit_n           = 50e-6;
	double n_p_cgs          = 1e17;
	double m_ion_amu        = 1.00794;
	double sz               = 30e-6;
	double sdelta           = 0.04;
	double t_tot            = 1.58631e-12*10;
	int n_steps             = 1;
	double dt               = t_tot/n_steps;
	std::string filename    = "output.h5";
	pushmethod_t pushmethod = ionsim::PUSH_SIMPLE;
	long n_field_x          = 21;
	long n_field_y          = 21;
	long n_field_z          = 11;

	const SimParams simparams(
		E,
		dt,
		emit_n,
		length,
		m_ion_amu,
		n_p_cgs,
		q_tot,
		radius,
		sdelta,
		sz,
		t_tot,
		n_steps,
		pushmethod,
		n_e,
		n_field_x,
		n_field_y,
		n_field_z,
		n_ions,
		filename
		);

	Field_Data *field;
	Field_Data *field_interp;
	Field_Comm fieldcomm;

	simparams.bcast_send();

	for (int step=0; step < n_steps; step++)
	{
		printf("Step: %d\n", step);
		field = new Field_Data(simparams);
		ionsim::loop_get_fields(fieldcomm, *field);

		/*
		switch (pushmethod)
		{
			case ionsim::PUSH_SIMPLE:
				ionsim::sendloop(ionsim::LOOP_PUSH_IONS);
				break;
			case ionsim::PUSH_FIELD:
				ionsim::loop_push_ions(*field);
				break;
		}
		*/

		/* (*field).dump_serial(simparams.filename, step); */

		/* field_interp = new Field(1001, 1001, simparams.radius, simparams.radius); */
		/* printf("Interp field init'ed\n"); */
		/* (*field).get_interp(*field_interp); */
		/* (*field_interp).dump_serial(simparams.filename, step); */
		/* delete field_interp; */

		delete field;

		/* ionsim::sendloop(ionsim::LOOP_DUMP_IONS, step); */
		/* ionsim::sendloop(ionsim::LOOP_DUMP_E, step); */
	}

	ionsim::sendloop(ionsim::LOOP_KILL);

	return 0;
}
