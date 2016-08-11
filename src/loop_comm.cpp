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
	MPI_Group world_group, slave_group;
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

int LoopComm::instruct(int *buf) const
{
	MPI_Bcast(buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
	return 0;
}

int LoopComm::instruct(const int buf) const
{
	int temp = buf;
	MPI_Bcast(&temp, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
	return 0;
}

int LoopComm::send_slaves(const int buf) const
{
	if (id != 0) return -1;

	int temp = buf;

	for (int i=1; i < p; i++)
	{
		MPI_Send(&temp, 1, MPI_INT, i, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);
	}

	return 0;
}

int LoopComm::recv_master(int *buf) const
{
	if (id == 0) return -1;

	MPI_Recv(buf, 1, MPI_INT, MASTER_RANK, TAG_LOOP_MESSAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	return 0;
}
