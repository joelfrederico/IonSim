#include "mpi.h"

int master(int &p)
{
	int root = 0;
	// ==============================
	// Starting gun
	// ==============================
	MPI::Status status;

	printf("I am the MASTER!\n");
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
	int n_e          = 1e6;
	int n_ion        = 1e6;
	double q_tot     = 2e10;
	double radius    = 100e-6;
	double length    = 100e-6;
	double E         = 20.35;
	double emit_n    = 50e-6;
	double n_p_cgs   = 1e18;
	double m_ion_amu = 1.00794;
	double sz        = 100e-6;
	double sdelta    = 0.04;

	MPI::COMM_WORLD.Bcast(&n_e       , 1 , MPI::LONG   , root);
	MPI::COMM_WORLD.Bcast(&n_ion     , 1 , MPI::LONG   , root);
	MPI::COMM_WORLD.Bcast(&radius    , 1 , MPI::DOUBLE , root);
	MPI::COMM_WORLD.Bcast(&length    , 1 , MPI::DOUBLE , root);
	MPI::COMM_WORLD.Bcast(&E         , 1 , MPI::DOUBLE , root);
	MPI::COMM_WORLD.Bcast(&emit_n    , 1 , MPI::DOUBLE , root);
	MPI::COMM_WORLD.Bcast(&n_p_cgs   , 1 , MPI::DOUBLE , root);
	MPI::COMM_WORLD.Bcast(&m_ion_amu , 1 , MPI::DOUBLE , root);
	MPI::COMM_WORLD.Bcast(&sz        , 1 , MPI::DOUBLE , root);
	MPI::COMM_WORLD.Bcast(&sdelta    , 1 , MPI::DOUBLE , root);
	
	// // ==============================
	// // Loop over sim steps
	// // ==============================
	// for (int i=0; i < max_sim_steps; i++)
	// {
	// 	// ==============================
	// 	// Loop over ion steps
	// 	// ==============================
	// 	for (int j=0; j < max_ion_steps; j++)
	// 	{
	// 		/* printf("Outer:\t%d; Inner:\t%d\n", i, j); */
	// 	}
	// }
	return 0;
}
