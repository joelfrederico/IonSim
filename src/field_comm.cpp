#include "field_comm.h"
#include "field_data.h"
#include "support_func.h"
#include <mpi.h>

Field_Comm::Field_Comm()
{
	p  = MPI::COMM_WORLD.Get_size();
	my_id = MPI::COMM_WORLD.Get_rank();
}

int Field_Comm::recv_field_others(Field_Data &field_recv)
{
	double *xbuf;
	double *ybuf;
	double *zbuf;

	xbuf = new double[field_recv.n_pts];
	ybuf = new double[field_recv.n_pts];
	zbuf = new double[field_recv.n_pts];

	for (int id=0; id < p; id++)
	{
		if (id != my_id) {
			std::cout << "Receiving " << field_recv.n_pts << " from: " << id << std::endl;
			MPI::COMM_WORLD.Recv(xbuf, field_recv.n_pts, MPI::DOUBLE, id, ionsim::TAG_FIELD);
			MPI::COMM_WORLD.Recv(ybuf, field_recv.n_pts, MPI::DOUBLE, id, ionsim::TAG_FIELD);
			MPI::COMM_WORLD.Recv(zbuf, field_recv.n_pts, MPI::DOUBLE, id, ionsim::TAG_FIELD);
			std::cout << "Finished receiving from: " << id << std::endl;

			for (int ind=0; ind < field_recv.n_pts; ind++)
			{
				field_recv.Ex_ind(ind) += xbuf[ind];
				field_recv.Ey_ind(ind) += ybuf[ind];
				field_recv.Ez_ind(ind) += zbuf[ind];
			}
		}
	}

	delete [] xbuf;
	delete [] ybuf;
	delete [] zbuf;

	return 0;
}

int Field_Comm::recv_field(Field_Data &field_recv, int sender_id)
{
	double *xbuf;
	double *ybuf;
	double *zbuf;

	xbuf = new double[field_recv.n_pts];
	ybuf = new double[field_recv.n_pts];
	zbuf = new double[field_recv.n_pts];

	MPI::COMM_WORLD.Recv(xbuf, field_recv.n_pts, MPI::DOUBLE, sender_id, ionsim::TAG_FIELD);
	MPI::COMM_WORLD.Recv(ybuf, field_recv.n_pts, MPI::DOUBLE, sender_id, ionsim::TAG_FIELD);
	MPI::COMM_WORLD.Recv(zbuf, field_recv.n_pts, MPI::DOUBLE, sender_id, ionsim::TAG_FIELD);

	for (int j=0; j < field_recv.n_pts; j++)
	{
		field_recv.Ex_ind(j) += xbuf[j];
		field_recv.Ey_ind(j) += ybuf[j];
		field_recv.Ez_ind(j) += zbuf[j];
	}

	delete [] xbuf;
	delete [] ybuf;
	delete [] zbuf;

	return 0;
}

int Field_Comm::send_field(Field_Data &field_send, int dest_id)
{
	std::cout << "Sending " << field_send.n_pts << " to: " << dest_id << std::endl;
	MPI::COMM_WORLD.Send(field_send.x_data, field_send.n_pts, MPI::DOUBLE, dest_id, ionsim::TAG_FIELD);
	MPI::COMM_WORLD.Send(field_send.y_data, field_send.n_pts, MPI::DOUBLE, dest_id, ionsim::TAG_FIELD);
	MPI::COMM_WORLD.Send(field_send.z_data, field_send.n_pts, MPI::DOUBLE, dest_id, ionsim::TAG_FIELD);
	std::cout << "Finished sending to: " << dest_id << std::endl;

	return 0;
}

