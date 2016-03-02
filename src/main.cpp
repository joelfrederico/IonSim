#include "main.h"
#include "config.h"

#include "master.h"
#include "slave.h"

#include "mpi.h"

const int INT_TAG = 1;

int main(int argc, char **argv)
{
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

	// ==============================
	// Create slave comm
	// ==============================
	int p_slave = p-1;
	int *ranks;
	MPI::Intracomm slave_comm_id;
	MPI::Group world_group_id = MPI::COMM_WORLD.Get_group();
	MPI::Group slave_group_id;

	ranks = new int[p_slave];
	for (int i=0; i < p; i++)
	{
		ranks[i] = i+1;
	}

	slave_group_id = world_group_id.Incl(p_slave, ranks);
	delete[] ranks;
	slave_comm_id = MPI::COMM_WORLD.Create(slave_group_id);

	if (id == 0)
	{
		master(p, id);
	}
	else
	{
		slave(p, id, slave_comm_id);
		printf("Process %d finished\n", id);
	}

	MPI::Finalize();

	return 0;
}
