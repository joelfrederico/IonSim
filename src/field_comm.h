#ifndef __FIELD_COMM_H_INCLUDED__
#define __FIELD_COMM_H_INCLUDED__

#include "field_data.h"
#include <mpi.h>

class Field_Comm
{
	private:
		int p;
		int my_id;

	public:
		Field_Comm();

		int recv_field_others_add(Field_Data &field_recv);
		int recv_field_copy(Field_Data &field_recv, int sender_id);
		int send_field(Field_Data &field_send, int dest_id);
		int sum_fields(Field_Data &field, MPI_Comm *comm_id);
};

#endif
