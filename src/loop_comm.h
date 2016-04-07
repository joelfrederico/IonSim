#ifndef __LOOP_COMM_H_INCLUDED__
#define __LOOP_COMM_H_INCLUDED__

#include <mpi.h>

class LoopComm
{
	private:
		MPI_Group world_group, slave_group;
		MPI_Comm _slave_create();

	public:
		const MPI_Comm slave_comm;
		const int p;
		const int id;

		LoopComm();
};

#endif
