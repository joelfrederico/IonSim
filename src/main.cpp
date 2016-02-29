#include "config.h"
#include "main.h"
#include "consts.h"

#include <stdio.h>
#include "mpi.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "hdf5.h"

#include "classes.h"

using namespace std;

int main(int argc, char **argv)
{
	// ==============================
	// Initialize MPI
	// ==============================
	int id;
	int p;
	MPI::Init ( argc, argv );

	// ==============================
	// Get the number of processes.
	// ==============================
	p = MPI::COMM_WORLD.Get_size ( );

	// ==============================
	// Get the individual process ID.
	// ==============================
	id = MPI::COMM_WORLD.Get_rank ( );

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

void master(int p)
{
	printf("I am the MASTER!\n");
	for (int i=1; i < p; i++)
	{
		MPI::COMM_WORLD.Send(&i, 1, MPI::INT, i, i*2);
	}
	// ==============================
	// Set up sim
	// ==============================
	int max_sim_steps = 100;
	int max_ion_steps = 100;
	
	// ==============================
	// Generate beam
	// ==============================
	
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
}

void slave(int id)
{
	int whoami = 0;
	int tag=1;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(&whoami, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
	tag = status.Get_tag();
	printf("Process %d got tag %d and says: ***DUM DUM DUM DUM***\n", whoami, tag);
}
