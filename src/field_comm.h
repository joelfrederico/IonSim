#ifndef __FIELD_COMM_H_INCLUDED__
#define __FIELD_COMM_H_INCLUDED__

#include "field_data.h"

class Field_Comm
{
	private:
		int p;
		int my_id;

	public:
		Field_Comm();

		int recv_field_others(Field_Data &field_recv);
		int recv_field(Field_Data &field_recv, int sender_id);
		int send_field(Field_Data &field_send, int dest_id);
};

#endif
