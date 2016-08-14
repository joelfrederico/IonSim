#ifndef __FFTW_CLASSES_H_INCLUDED__
#define __FFTW_CLASSES_H_INCLUDED__

#include "scalar_data.h"
#include <fftw3-mpi.h>

// ==================================
// Functions for sending data
// ==================================
int psifftw_base(SimParams simparams, LoopComm loopcomm);
int fftwl_recv_local_size(long long &local_n0, long long &local_0_start, const int id);

template<typename T>
int MPI_Recv_complex(int id, std::vector<std::complex<T>> &buf)
{
	int count;
	long long local_0_start_ll;
	long ind;
	MPI_Status status;
	MPI_Datatype mpitype;
	std::vector<T> dbuf;
	
	MPI_Recv(&local_0_start_ll, 1, MPI_LONG_LONG, 1, TAG_COMPLEX_VEC_START, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	MPI_Probe(id, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	if (status.MPI_TAG == TAG_LDOUBLE_COMPLEX_VEC)
	{
		JTF_PRINT(Receiving long double);
		mpitype = MPI_LONG_DOUBLE;
	} else if (status.MPI_TAG == TAG_DOUBLE_COMPLEX_VEC) {
		mpitype = MPI_DOUBLE;
	} else if (status.MPI_TAG == TAG_FLOAT_COMPLEX_VEC) {
		mpitype = MPI_FLOAT;
	}

	MPI_Get_count(&status, mpitype, &count);
	JTF_PRINT_NOEND(Receiving count: ) << count << std::endl;
	dbuf.resize(count);
	buf.resize(count/2);

	MPI_Recv(dbuf.data(), count, mpitype, id, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	for (long i=0; i<count/2; i++)
	{
		ind = i*2;
		buf[i + local_0_start_ll].real(dbuf[ind]);
		buf[i].imag(dbuf[ind+1]);
	}

	return 0;
}

// ==================================
// Classes: receiving and doing FFT
// ==================================
class Fftw_base
{
	private:

	public:
		Fftw_base(const ptrdiff_t N0, const ptrdiff_t N1, LoopComm loopcomm);
};

class Fftw_2d_dft
{
	private:

	public:
		Fftw_2d_dft();

		int receive_data();
		int send_data();
		int execute();
};


class Fftw_2d_idft
{
	private:

	public:

};

#endif
