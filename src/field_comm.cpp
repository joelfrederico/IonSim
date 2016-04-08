#include "field_comm.h"
#include "field_data.h"
#include "support_func.h"
#include <mpi.h>
#include "loop_comm.h"
#include "consts.h"

Field_Comm::Field_Comm()
{
	LoopComm loopcomm;
	p     = loopcomm.p;
	my_id = loopcomm.id;
}

int Field_Comm::recv_field_others_add(Field_Data &field_recv)
{
	MPI_Status status;
	double *xbuf;
	double *ybuf;
	double *zbuf;

	xbuf = new double[field_recv.n_pts];
	ybuf = new double[field_recv.n_pts];
	zbuf = new double[field_recv.n_pts];

	for (int id=0; id < p; id++)
	{
		if (id != my_id)
		{
			/* std::cout << "Receiving " << field_recv.n_pts << " from: " << id << std::endl; */
			MPI_Recv(xbuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			MPI_Recv(ybuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			MPI_Recv(zbuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			/* std::cout << "Finished receiving from: " << id << std::endl; */

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

int Field_Comm::recv_field_copy(Field_Data &field_recv, int sender_id)
{
	MPI_Status status;
	double *xbuf;
	double *ybuf;
	double *zbuf;

	xbuf = new double[field_recv.n_pts];
	ybuf = new double[field_recv.n_pts];
	zbuf = new double[field_recv.n_pts];

	MPI_Recv(xbuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);
	MPI_Recv(ybuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);
	MPI_Recv(zbuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);

	for (int j=0; j < field_recv.n_pts; j++)
	{
		field_recv.Ex_ind(j) = xbuf[j];
		field_recv.Ey_ind(j) = ybuf[j];
		field_recv.Ez_ind(j) = zbuf[j];
	}

	delete [] xbuf;
	delete [] ybuf;
	delete [] zbuf;

	return 0;
}

int Field_Comm::send_field(Field_Data &field_send, int dest_id)
{
	/* std::cout << "Sending " << field_send.n_pts << " to: " << dest_id << std::endl; */
	MPI_Send(field_send.x_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	MPI_Send(field_send.y_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	MPI_Send(field_send.z_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	/* std::cout << "Finished sending to: " << dest_id << std::endl; */

	return 0;
}

int Field_Comm::sum_fields(Field_Data &field, MPI_Comm *comm_id)
{
	double *xbuf;
	double *ybuf;
	double *zbuf;

	xbuf = new double[field.n_pts];
	ybuf = new double[field.n_pts];
	zbuf = new double[field.n_pts];

	Field_Data temp(field.x_pts, field.y_pts, field.z_pts, field.x_edge_mag, field.y_edge_mag, field.z_edge_mag);;

	for (int i=0; i < p; i++)
	{
		if (i == my_id)
		{
			MPI_Bcast(field.x_data, field.n_pts, MPI_DOUBLE, 0, *comm_id);
		} else {
			MPI_Bcast(&xbuf, field.n_pts, MPI_DOUBLE, 0, *comm_id);
		}
	}
	return 0;
}
