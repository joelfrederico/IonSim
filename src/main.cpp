#include "config.h"
#include "classes.h"
#include "consts.h"
#include "main.h"

#include <stdio.h>
#include "mpi.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "hdf5.h"


const int INT_TAG = 1;

int main(int argc, char **argv)
{
	float E = 20.35;
	Emit emit(50e-6, E, true);
	Plasma plas(1e18, 1.00794);
	Match mat(plas, E, emit);

	printf("n_p: %.3e, m: %.3e, w_p: %.3e, k_ion(20.35): %.3e\n", plas.n_p(), plas.m(), plas.w_p(), plas.k_ion(20.35));
	printf("Match- beta: %.6e\n", mat.beta());

	// ==============================
	// Initialize MPI
	// ==============================
	int id;
	int p;
	MPI::Init(argc, argv);

	// ==============================
	// Get the number of processes.
	// ==============================
	p = MPI::COMM_WORLD.Get_size();

	// ==============================
	// Get the individual process ID.
	// ==============================
	id = MPI::COMM_WORLD.Get_rank();

	if (id == 0)
	{
		master(p);
	}
	else
	{
		slave(id);
	}

	MPI::Finalize();

	return 0;
}

int master(int p)
{
	int output;
	int tag;


	MPI::Status status;

	printf("I am the MASTER!\n");
	for (int id=1; id < p; id++)
	{
		tag = id*2;
		MPI::COMM_WORLD.Send(&id, 1, MPI::INT, id, id*2);
		MPI::COMM_WORLD.Recv(&output, 1, MPI::INT, id, MPI::ANY_TAG, status);
		printf("The Master received %d\n", output);
	}
	// ==============================
	// Set up sim
	// ==============================
	int max_sim_steps = 100;
	int max_ion_steps = 100;
	
	// ==============================
	// Generate beam
	// ==============================
	Beam bea (1, 0);
	
	// ==============================
	// Generate ions
	// ==============================
	
	// ==============================
	// Loop over sim steps
	// ==============================
	for (int i=0; i < max_sim_steps; i++)
	{
		// ==============================
		// Loop over ion steps
		// ==============================
		for (int j=0; j < max_ion_steps; j++)
		{
			/* printf("Outer:\t%d; Inner:\t%d\n", i, j); */
		}
	}
	return 0;
}

int slave(int id)
{
	int whoami = 0;
	int tag    = 1;
	int output = 3;

	MPI::Status status;

	MPI::COMM_WORLD.Recv(&whoami, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
	tag = status.Get_tag();
	printf("Process %d got tag %d and says: ***DUM DUM DUM DUM***\n", whoami, tag);

	output *= tag;
	
	MPI::COMM_WORLD.Send(&output, 1, MPI::INT, 0, INT_TAG); 
	return 0;
}
