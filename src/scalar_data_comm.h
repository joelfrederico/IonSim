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
		int recv_scalar_others_add(ScalarData<T> &scalar_recv) const;

		template<typename T>
		int recv_scalar_copy(ScalarData<T> &scalar_recv, const int sender_id) const ;

		template<typename T>
		int send_scalar(const ScalarData<T> &scalar_send, const int dest_id) const ;

		template<typename T>
		int bcast_scalar(const ScalarData<T> &scalar_bcast) const;
};

#endif
