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

		template<typename T>
		int recv_scalar_others_add(ScalarData<T> &scalar_recv);

		template<typename T>
		int recv_scalar_copy(ScalarData<T> &scalar_recv, int sender_id);

		template<typename T>
		int send_scalar(ScalarData<T> &scalar_send, int dest_id);

		template<typename T>
		int bcast_scalar(ScalarData<T> &scalar_bcast);
};

#endif
