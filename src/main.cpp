#include "config.h"
#include "main.h"

#include "master.h"
#include "slave.h"

#include <mpi.h>
#include <unistd.h>

const int INT_TAG = 1;

int main(int argc, char **argv)
{
	int id, p, c;
	int p_slave;
	int *ranks;
	MPI::Group world_group_id, slave_group_id;
	MPI::Intracomm slave_comm_id;
	bool verbose;

	verbose = false;

	/* while ((c = getopt (argc, argv, "v")) != -1) */
	/* { */
	/* 	switch (c) */
	/* 		{ */
	/* 		case 'v': */
	/* 			verbose = true; */
	/* 			break; */
	/* 		default: */
	/* 			abort (); */
	/* 		} */
	/* } */

	// ==============================
	// Initialize MPI
	// ==============================
	MPI::Init(argc, argv);

	// ==============================
	// Get the number of processes.
	// ==============================
	p = MPI::COMM_WORLD.Get_size();

	// ==============================
	// Get the individual process ID.
	// ==============================
	id = MPI::COMM_WORLD.Get_rank();

	if (verbose && id==0) printf("%d %d %d\n", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION);
	// ==============================
	// Create slave comm
	// ==============================
	p_slave = p-1;
	world_group_id = MPI::COMM_WORLD.Get_group();

	ranks = new int[p_slave];
	for (int i=0; i < p_slave; i++)
	{
		ranks[i] = i+1;
	}

	slave_group_id = world_group_id.Incl(p_slave, ranks);
	delete[] ranks;
	slave_comm_id = MPI::COMM_WORLD.Create(slave_group_id);

	// ==============================
	// Simulate things
	// ==============================
	if (id == 0) {
		master(p, verbose);
	} else {
		slave(p, id, slave_comm_id, verbose);
	}

	MPI::Finalize();

	if (verbose)
	{
		if (id == 0)
		{
			printf("Master finished\n");
		}
		else
		{
			printf("Slave %d finished\n", id);
		}
	}


	return 0;
}
