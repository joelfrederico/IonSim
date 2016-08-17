#include <complex>
#include "fftw_classes.h"
#include "loop_comm.h"
#include "scalar_data.h"
#include <fftw3-mpi.h>
#include <mpi.h>
#include "support_func.h"
#include "mpi_vec.h"
#include <math.h>

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

double f(unsigned long long i, unsigned long long j)
{
	long double valx = j;
	long double valy = i;
	valx = (valx-(256-1)/2);
	valy = (valy-(256-1)/2);

	if (((i == 127) || (i == 128)) && ((j==127) || (j==128)))
	{
		valx = 0.25;
	} else {
		valx = 0;
	}
	return valx;
}

int psifftw_base(LoopComm loopcomm, ptrdiff_t N0, ptrdiff_t N1)
{
	// ==================================
	// Initialize Variables
	// ==================================
	long double *r_buf;
	fftwl_complex *c_buf;
	ptrdiff_t alloc_local, local_n0, local_0_start, i_nonlocal;
	fftwl_plan planforward, planreverse;
	std::vector<std::complex<long double>> cdata;
	long ind;
	long double delx, dely;
	delx = dely = 0.05;

	// ==================================
	// Load parameters
	// ==================================
	N0 = 256;
	N1 = 256;

	// ==================================
	// Set up FFT
	// ==================================
	alloc_local = fftwl_mpi_local_size_2d(N0, int(N1/2+1), loopcomm.slave_comm, &local_n0, &local_0_start);

	c_buf = fftwl_alloc_complex(alloc_local);
	r_buf = fftwl_alloc_real(2*alloc_local);

	planforward = fftwl_mpi_plan_dft_r2c_2d(N0, N1, r_buf, c_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE);

	// Copy in data
	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		i_nonlocal = i+local_0_start;
		for (ptrdiff_t j=0; j<N1; j++)
		{
			ind = i*(N1+2)+j;
			r_buf[ind] = f(i_nonlocal, j);
		}
	}

	// ==================================
	// Perform FFT
	// ==================================
	fftwl_execute(planforward);

	// ==================================
	// Copy FFT results
	// ==================================
	cdata = MPI_convert_complex_vec(c_buf, local_n0, N1);

	// ==================================
	// Divide by k^2
	// ==================================
	MPI_Complex_div_k2(cdata, delx, dely, local_n0, local_0_start, N0, N1);

	// ==================================
	// Send divided results to master
	// ==================================
	MPI_Send_Scalar_Complex(cdata, local_n0, local_0_start, N0, N1);

	// ==================================
	// Set up IFFT
	// ==================================
	planreverse = fftwl_mpi_plan_dft_c2r_2d(N0, N1, c_buf, r_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE);

	// Copy in data
	MPI_convert_complex_buf(cdata, c_buf, local_n0, local_0_start, N1);

	// ==================================
	// Perform IFFT
	// ==================================
	fftwl_execute(planreverse);

	// ==================================
	// Copy IFFT results
	// ==================================
	auto rdata = MPI_convert_real_vec(r_buf, local_n0, N0, N1);

	// ==================================
	// Send results to master
	// ==================================
	MPI_Send_Scalar_Real(rdata, local_n0, local_0_start, N0, N1);

	// ==================================
	// Deallocate
	// ==================================
	fftwl_free(r_buf);
	fftwl_free(c_buf);
	fftwl_destroy_plan(planforward);
	fftwl_destroy_plan(planreverse);

	return 0;
}
