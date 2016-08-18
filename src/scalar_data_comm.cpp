#include "scalar_data_comm.h"
#include "scalar_data.h"
#include "support_func.h"
#include <mpi.h>
#include "loop_comm.h"
#include "consts.h"
#include <vector>

ScalarData_Comm::ScalarData_Comm()
{
	LoopComm loopcomm;
	p     = loopcomm.p;
	my_id = loopcomm.id;
}

template<typename T>
int ScalarData_Comm::recv_scalar_others_add(ScalarData<T> &scalar_recv)
{
	MPI_Status status;
	int count;
	int n_pts = scalar_recv.n_pts();

	std::vector<long double> buf;
	buf.resize(n_pts);

	for (int id=0; id < p; id++)
	{
		if (id != my_id)
		{
			/* std::cout << "Receiving " << n_pts << " from: " << id << std::endl; */

			MPI_Probe(id, TAG_FIELD, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_LONG_DOUBLE, &count);
			if (count != n_pts) throw std::runtime_error("Receiving size does not match sending size.");

			MPI_Recv(buf.data(), n_pts, MPI_LONG_DOUBLE, id, TAG_FIELD, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			/* std::cout << "Finished receiving from: " << id << std::endl; */

			for (typename decltype(scalar_recv.vdata())::size_type ind=0; ind < n_pts; ind++)
			{
				scalar_recv.ind(ind) += buf[ind];
			}
		}
	}

	return 0;
}

template<typename T>
int ScalarData_Comm::recv_scalar_copy(ScalarData<T> &scalar_recv, int sender_id)
{
	auto n_pts = scalar_recv.n_pts();

	MPI_Recv(scalar_recv.data.data(), n_pts, MPI_DOUBLE, sender_id, TAG_FIELD, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	return 0;
}

template<typename T>
int ScalarData_Comm::bcast_scalar(ScalarData<T> &scalar_bcast)
{
	MPI_Datatype mpi_type;
	unsigned long n_pts = scalar_bcast.n_pts();
	LoopComm loopcomm;

	if (typeid(T) == typeid(long double))
	{
		mpi_type = MPI_LONG_DOUBLE;
	} else {
		throw std::runtime_error("Type not supported.");
	}

	MPI_Bcast(&n_pts, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if (n_pts != scalar_bcast.n_pts()) throw std::runtime_error("Incompatible number of points.");

	MPI_Bcast(scalar_bcast.vdata().data(), n_pts, mpi_type, 0, MPI_COMM_WORLD);

	return 0;
}

template<typename T>
int ScalarData_Comm::send_scalar(ScalarData<T> &scalar_send, int dest_id)
{
	int n_pts = scalar_send.n_pts();
	auto vdata = scalar_send.vdata();

	/* std::cout << "Sending " << n_pts << " to: " << dest_id << std::endl; */
	if (vdata.size() != n_pts) throw std::runtime_error("Size to send does not match vector.");
	MPI_Send(vdata.data(), n_pts, MPI_LONG_DOUBLE, dest_id, TAG_FIELD, MPI_COMM_WORLD);
	/* std::cout << "Finished sending to: " << dest_id << std::endl; */

	return 0;
}

template int ScalarData_Comm::send_scalar(ScalarData<ldouble> &scalar_send, int dest_id);
template int ScalarData_Comm::recv_scalar_others_add(ScalarData<ldouble> &scalar_recv);
template int ScalarData_Comm::bcast_scalar(ScalarData<long double> &scalar_bcast);
