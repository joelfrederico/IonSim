#ifndef __SCALAR_DATA_COMM_H_INCLUDED__
#define __SCALAR_DATA_COMM_H_INCLUDED__

#include "scalar_data.h"
#include <mpi.h>

class ScalarData_Comm
{
	private:
		int p;
		int my_id;

	public:
		ScalarData_Comm();

		int recv_scalar_others_add(ScalarData &scalar_recv);
		int recv_scalar_copy(ScalarData &scalar_recv, int sender_id);
		int send_scalar(ScalarData &scalar_send, int dest_id);
};

#endif
