#include "loop_comm.h"
#include <mpi.h>

int _p()
{
	int out;
	MPI_Comm_size(MPI_COMM_WORLD, &out);
	return out;
}

int _id()
{
	int out;
	MPI_Comm_rank(MPI_COMM_WORLD, &out);
	return out;
}

MPI_Comm LoopComm::_slave_create()
{
	int p = _p();
	int id = _id();
	MPI_Comm out;
	int status;

	// ==============================
	// Create slave comm
	// ==============================
	int p_slave;
	int *ranks;

	p_slave = p-1;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	ranks = new int[p_slave];
	for (int i=0; i < p_slave; i++)
	{
		ranks[i] = i+1;
	}

	MPI_Group_incl(world_group, p_slave, ranks, &slave_group);
	delete[] ranks;

	status = MPI_Comm_create_group(MPI_COMM_WORLD, slave_group, 0, &out);

	return out;
}

LoopComm::LoopComm() :
	p          ( _p()            ) ,
	id         ( _id()           ) ,
	slave_comm ( _slave_create() )
{
}
