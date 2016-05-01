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
	double *Exbuf;
	double *Eybuf;
	double *Ezbuf;

	double *Bxbuf;
	double *Bybuf;
	double *Bzbuf;

	Exbuf = new double[field_recv.n_pts];
	Eybuf = new double[field_recv.n_pts];
	Ezbuf = new double[field_recv.n_pts];

	Bxbuf = new double[field_recv.n_pts];
	Bybuf = new double[field_recv.n_pts];
	Bzbuf = new double[field_recv.n_pts];

	for (int id=0; id < p; id++)
	{
		if (id != my_id)
		{
			/* std::cout << "Receiving " << field_recv.n_pts << " from: " << id << std::endl; */
			MPI_Recv(Exbuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			MPI_Recv(Eybuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			MPI_Recv(Ezbuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);

			MPI_Recv(Bxbuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			MPI_Recv(Bybuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			MPI_Recv(Bzbuf, field_recv.n_pts, MPI_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, &status);
			/* std::cout << "Finished receiving from: " << id << std::endl; */

			for (int ind=0; ind < field_recv.n_pts; ind++)
			{
				field_recv.Ex_ind(ind) += Exbuf[ind];
				field_recv.Ey_ind(ind) += Eybuf[ind];
				field_recv.Ez_ind(ind) += Ezbuf[ind];

				field_recv.Ex_ind(ind) += Bxbuf[ind];
				field_recv.Ey_ind(ind) += Bybuf[ind];
				field_recv.Ez_ind(ind) += Bzbuf[ind];
			}
		}
	}

	delete [] Exbuf;
	delete [] Eybuf;
	delete [] Ezbuf;

	delete [] Bxbuf;
	delete [] Bybuf;
	delete [] Bzbuf;

	return 0;
}

int Field_Comm::recv_field_copy(Field_Data &field_recv, int sender_id)
{
	MPI_Status status;
	double *Exbuf;
	double *Eybuf;
	double *Ezbuf;

	double *Bxbuf;
	double *Bybuf;
	double *Bzbuf;

	Exbuf = new double[field_recv.n_pts];
	Eybuf = new double[field_recv.n_pts];
	Ezbuf = new double[field_recv.n_pts];

	Bxbuf = new double[field_recv.n_pts];
	Bybuf = new double[field_recv.n_pts];
	Bzbuf = new double[field_recv.n_pts];

	MPI_Recv(Exbuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);
	MPI_Recv(Eybuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);
	MPI_Recv(Ezbuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);

	MPI_Recv(Bxbuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);
	MPI_Recv(Bybuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);
	MPI_Recv(Bzbuf, field_recv.n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, &status);

	for (int j=0; j < field_recv.n_pts; j++)
	{
		field_recv.Ex_ind(j) = Exbuf[j];
		field_recv.Ey_ind(j) = Eybuf[j];
		field_recv.Ez_ind(j) = Ezbuf[j];

		field_recv.Ex_ind(j) = Bxbuf[j];
		field_recv.Ey_ind(j) = Bybuf[j];
		field_recv.Ez_ind(j) = Bzbuf[j];
	}

	delete [] Exbuf;
	delete [] Eybuf;
	delete [] Ezbuf;

	delete [] Bxbuf;
	delete [] Bybuf;
	delete [] Bzbuf;

	return 0;
}

int Field_Comm::send_field(Field_Data &field_send, int dest_id)
{
	/* std::cout << "Sending " << field_send.n_pts << " to: " << dest_id << std::endl; */
	MPI_Send(field_send.Ex_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	MPI_Send(field_send.Ey_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	MPI_Send(field_send.Ez_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);

	MPI_Send(field_send.Bx_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	MPI_Send(field_send.By_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	MPI_Send(field_send.Bz_data, field_send.n_pts, MPI_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	/* std::cout << "Finished sending to: " << dest_id << std::endl; */

	return 0;
}

int Field_Comm::sum_fields(Field_Data &field, MPI_Comm *comm_id)
{
	double *Exbuf;
	double *Eybuf;
	double *Ezbuf;

	double *Bxbuf;
	double *Bybuf;
	double *Bzbuf;

	Exbuf = new double[field.n_pts];
	Eybuf = new double[field.n_pts];
	Ezbuf = new double[field.n_pts];

	Bxbuf = new double[field.n_pts];
	Bybuf = new double[field.n_pts];
	Bzbuf = new double[field.n_pts];

	Field_Data temp(field.x_pts, field.y_pts, field.z_pts, field.x_edge_mag, field.y_edge_mag, field.z_edge_mag);;

	for (int i=0; i < p; i++)
	{
		if (i == my_id)
		{
			MPI_Bcast(field.Ex_data, field.n_pts, MPI_DOUBLE, 0, *comm_id);
		} else {
			MPI_Bcast(&Exbuf, field.n_pts, MPI_DOUBLE, 0, *comm_id);
		}
	}
	return 0;

	delete [] Exbuf;
	delete [] Eybuf;
	delete [] Ezbuf;

	delete [] Bxbuf;
	delete [] Bybuf;
	delete [] Bzbuf;
}
