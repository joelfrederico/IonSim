#include <complex>
#include "fftw_classes.h"
#include "loop_comm.h"
#include "scalar_data.h"
#include <fftw3-mpi.h>
#include <mpi.h>

int fftwl_send_local_size(ptrdiff_t &local_n0, ptrdiff_t &local_0_start)
{
	// ==================================
	// Tell Master local data size
	// ==================================
	MPI_Send(&local_n0, 1, MPI_LONG_LONG_INT, 0, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);
	MPI_Send(&local_0_start, 1, MPI_LONG_LONG_INT, 0, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);

	return 0;
}

int fftwl_recv_local_size(long long &local_n0, long long &local_0_start, const int id)
{
	// Get array size, starting index
	MPI_Recv(&local_n0, 1, MPI_LONG_LONG_INT, id, TAG_LOOP_MESSAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&local_0_start, 1, MPI_LONG_LONG_INT, id, TAG_LOOP_MESSAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	return 0;
}

template<typename T>
int MPI_Send_complex(const std::vector<std::complex<T>> &buf)
{
	typename std::vector<T> dbuf;
	long ind;
	int count;
	MPI_Datatype mpitype;
	int tag;

	count = buf.size();
	dbuf.resize(2*count);
	for (long i=0; i<count; i++)
	{
		ind = i*2;
		dbuf[ind]   = buf[i].real();
		dbuf[ind+1] = buf[i].imag();
	}

	if (typeid(T) == typeid(long double))
	{
		mpitype = MPI_LONG_DOUBLE;
		tag = TAG_LDOUBLE_COMPLEX_VEC;
	} else if (typeid(T) == typeid(double)) {
		mpitype = MPI_DOUBLE;
		tag = TAG_DOUBLE_COMPLEX_VEC;
	} else if (typeid(T) == typeid(float)) {
		mpitype = MPI_FLOAT;
		tag = TAG_FLOAT_COMPLEX_VEC;
	} else {
		throw std::runtime_error("Not a valid complex type!");
	}

	MPI_Send(dbuf.data(), 2*count, mpitype, 0, tag, MPI_COMM_WORLD);

	return 0;
}


long c_ind(const long i, const long j, const ptrdiff_t N1)
{
	return i*(2*(N1/2+1)) + j;
}

long column_major(const long i, const long j, const ptrdiff_t N0)
{
	return i + N0*j;
}

long row_major(const long i, const long j, const ptrdiff_t N1)
{
	return i*N1 + j;
}

int psifftw_base(SimParams simparams, LoopComm loopcomm)
{
	// ==================================
	// Allocate variables
	// ==================================
	fftwl_plan planforward, planreverse;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	ldouble_vec rho_x;
	long long int rho_x_size;
	long long int i, j;
	ScalarData psi(simparams);
	long double f0, f1, kx, ky, k2;
	long double *real_out, *r_buf;
	fftwl_complex *c_buf, *rho_k;
	long ind, ind2;
	std::vector<std::complex<long double>> cdata;

	// ==================================
	// Size of FFT
	// ==================================
	ptrdiff_t N0 = simparams.n_field_x;
	ptrdiff_t N1 = simparams.n_field_y;

	// ==================================
	// Determine local size
	// ==================================
	alloc_local = fftwl_mpi_local_size_2d(N0, N1/2+1, loopcomm.slave_comm, &local_n0, &local_0_start);

	// ==================================
	// Receive local data
	// ==================================
	fftwl_send_local_size(local_n0, local_0_start);
	rho_x_size = local_n0*N1;
	rho_x.resize(rho_x_size);
	MPI_Recv(rho_x.data(), rho_x_size, MPI_LONG_DOUBLE, 0, TAG_LOOP_MESSAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// ==================================
	// Allocate complex data
	// LARGER THAN NECESSARY
	// ==================================
	rho_k = fftwl_alloc_complex(alloc_local);
	c_buf = fftwl_alloc_complex(alloc_local);
	r_buf = fftwl_alloc_real(2*alloc_local);
	real_out = new long double[local_n0*N1];

	// ==================================
	// Create and execute plan
	// ==================================
	planforward = fftwl_mpi_plan_dft_r2c_2d(N0, N1, r_buf, c_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);

	for (long i=0; i < local_n0; i++)
	{
		for (long j=0; j < N1; j++)
		{
			r_buf[i*(2*(N1/2+1))+j] = rho_x[i + N0*j];
		}
	}

	fftwl_execute(planforward);

	cdata.resize(N0*(N1/2+1));
	// ==================================
	// Divide by k^2
	// ==================================
	f0 = psi.dxdi * 2 * N0;
	f1 = psi.dydj * 2 * N1;
	for (long i=0; i < local_n0; i++)
	{
		if (i <= N0/2)
		{
			kx = i / f0;
		} else {
			kx = (i-N0) / f0;
		}
		for (long j=0; j< N1/2+1; j++)
		{
			ky = j/f1;
			k2 = kx*kx + ky*ky;
			ind = c_ind(i, j, N1);
			rho_k[ind][0] = c_buf[ind][0] / k2;
			rho_k[ind][1] = c_buf[ind][1] / k2;
			ind2 = row_major(i, j, N0);
			cdata[ind2].real(rho_k[ind][0]);
			cdata[ind2].imag(rho_k[ind][1]);
		}
	}

	MPI_Send_complex(cdata);

	// ==================================
	// Create and execute plan
	// ==================================
	JTF_PRINTVAL_NOEND(alloc_local) << ", id: " << loopcomm.id << std::endl;
	JTF_PRINTVAL_NOEND(r_buf)    << ", size: " << sizeof(long double)   << ", id: " << loopcomm.id << std::endl;
	JTF_PRINTVAL_NOEND(c_buf)    << ", size: " << sizeof(fftwl_complex) << ", id: " << loopcomm.id << std::endl;
	JTF_PRINTVAL_NOEND(real_out) << ", size: " << sizeof(long double)   << ", id: " << loopcomm.id << std::endl;
	JTF_PRINTVAL_NOEND(rho_k)    << ", size: " << sizeof(fftwl_complex) << ", id: " << loopcomm.id << std::endl;
	JTF_PRINT(Or die trying);
	fftwl_free(r_buf);
	fftwl_free(rho_k);
	fftwl_free(c_buf);
	fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);
	delete[] real_out;
	return 0;

	fftwl_free(c_buf);
	fftwl_free(r_buf);

	alloc_local = fftwl_mpi_local_size_2d(N0, N1/2+1, loopcomm.slave_comm, &local_n0, &local_0_start);
	c_buf = fftwl_alloc_complex(alloc_local);
	r_buf = fftwl_alloc_real(2*alloc_local);

	planreverse = fftwl_mpi_plan_dft_c2r_2d(N0, N1, c_buf, r_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
	JTF_PRINT(Or try dying);

	for (long i=0; i < N0; i++)
	{
		for (long j=0; j<N1/2+1; j++)
		{
			c_buf[c_ind(i, j, N1)][0] = rho_k[row_major(i, j, N1)][0];
			c_buf[c_ind(i, j, N1)][1] = rho_k[row_major(i, j, N1)][1];
		}
	}

	JTF_PRINT(Yaya);

	fftwl_execute(planreverse);

	JTF_PRINT(Nana);
	for (long i=0; i < local_n0; i++)
	{
		for (long j=0; j<N1; j++)
		{
			real_out[column_major(i, j, local_n0)] = r_buf[row_major(i, j, N1)];
		}
	}

	JTF_PRINT(Further);

	// ==================================
	// Send Data
	// ==================================
	fftwl_send_local_size(local_n0, local_0_start);
	MPI_Send(real_out, rho_x_size, MPI_LONG_DOUBLE, 0, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);

	fftwl_free(r_buf);
	fftwl_free(rho_k);
	delete[] real_out;

	// ==================================
	// Send wisdom back to master
	// ==================================
	fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);

	fftwl_destroy_plan(planforward);
	fftwl_destroy_plan(planreverse);

	return 0;
}
