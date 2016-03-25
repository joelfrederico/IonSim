#include "mpi.h"
#include "support_func.h"

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
	double radius           = 2.4276628847185805e-06*10;
	double length           = 100e-6;
	double E                = 20.35;
	double emit_n           = 50e-6;
	double n_p_cgs          = 1e17;
	double m_ion_amu        = 1.00794;
	double sz               = 30e-6;
	double sdelta           = 0.04;
	double t_tot            = 1.58631e-12*2;
	int n_steps             = 100;
	double dt               = t_tot/n_steps;
	std::string filename    = "output.h5";
	pushmethod_t pushmethod = ionsim::PUSH_FIELD;
	long n_field_x          = 101;
	long n_field_y          = 101;
	long n_field_z          = 101;

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

	Field *field;

	simparams.bcast_send();

	for (int step=0; step < n_steps; step++)
	{
		field = new Field(simparams);
		ionsim::loop_get_fields(*field);

		ionsim::loop_push_ions(*field);
		printf("Master field: %e\n", (*field).Ex_ind(0, 0));

		(*field).dump_serial(simparams.filename, step);

		delete field;

		ionsim::sendloop(ionsim::LOOP_DUMP_IONS, step);
		ionsim::sendloop(ionsim::LOOP_DUMP_E, step);
	}

	ionsim::sendloop(ionsim::LOOP_KILL);

	return 0;
}

int push_ions()
{
	Field *field;

	return 0;
}
