#include "mpi.h"
#include "support_func.h"

int master(int &p)
{
	// ==============================
	// Starting gun
	// ==============================
	MPI::Status status;

	printf("Master says: I am the MASTER!\n");
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
	
	long n_e             = 1e4;
	long n_ion           = 1e4;
	double q_tot         = 2e10;
	double radius        = 2.440175e-7*10;
	double length        = 100e-6;
	double E             = 20.35;
	double emit_n        = 50e-6;
	double n_p_cgs       = 1e17;
	double m_ion_amu     = 1.00794;
	double sz            = 30e-6;
	double sdelta        = 0.04;
	double t_tot         = 1.58631e-12*2;
	int n_steps          = 100;
	double dt            = t_tot/n_steps;
	std::string filename = "output.hdf5";
	int runge_kutta      = 0;
	long n_field_x       = 101;
	long n_field_y       = 101;
	long n_field_z       = 101;

	// Send numerical parameters
	MPI::COMM_WORLD.Bcast(&n_e         , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&n_ion       , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&q_tot       , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&radius      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&length      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&E           , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&emit_n      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&n_p_cgs     , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&m_ion_amu   , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&sz          , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&sdelta      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&t_tot       , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&n_steps     , 1 , MPI::INT    , 0);
	MPI::COMM_WORLD.Bcast(&dt          , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&runge_kutta , 1 , MPI::INT    , 0);
	MPI::COMM_WORLD.Bcast(&n_field_x   , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&n_field_y   , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&n_field_z   , 1 , MPI::LONG   , 0);

	// Send string
	char *cbuf;
	long cbuf_l;

	cbuf_l = filename.length()+1;
	MPI::COMM_WORLD.Bcast(&cbuf_l, 1, MPI::LONG, 0);
	cbuf = new char[cbuf_l];
	strcpy(cbuf, filename.c_str());
	MPI::COMM_WORLD.Bcast(cbuf, cbuf_l, MPI::CHAR, 0);
	delete [] cbuf;

	/* n_steps = 400; */
	for (int step=0; step < n_steps; step++)
	{
		ionsim::sendloop(&ionsim::LOOP_GET_EFIELD);
		ionsim::sendloop(&ionsim::LOOP_PUSH_IONS);
		ionsim::sendloop(&ionsim::LOOP_DUMP_IONS, step);
	}

	ionsim::sendloop(&ionsim::LOOP_KILL);

	return 0;
}
