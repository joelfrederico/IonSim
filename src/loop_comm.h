#ifndef __LOOP_COMM_H_INCLUDED__
#define __LOOP_COMM_H_INCLUDED__

#include <mpi.h>
#include "consts.h"

class LoopComm
{
	private:
		MPI_Comm _slave_create();

	public:
		const MPI_Comm slave_comm;
		const int p;
		const int id;

		LoopComm();

		int instruct(int *buf) const;
		int instruct(const int buf) const;

		int send_slaves(const int buf) const;
		int recv_master(int *buf) const;
};

#endif
