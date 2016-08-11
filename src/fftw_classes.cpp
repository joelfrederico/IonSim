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

int psifftw_base(SimParams simparams, LoopComm loopcomm)
{
	// ==================================
	// Allocate variables
	// ==================================
	fftwl_plan plan;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	ldouble_vec rho_x;
	long long int rho_x_size;
	fftwl_complex *rho_k;
	long long int i, j;
	ScalarData psi(simparams);
	long double f0, f1, kx, ky, k2;
	long double *real_out;

	// ==================================
	// Size of FFT
	// ==================================
	ptrdiff_t N0 = simparams.n_field_x;
	ptrdiff_t N1 = simparams.n_field_y;

	// ==================================
	// Determine local size
	// ==================================
	alloc_local = fftwl_mpi_local_size_2d(N0, N1, loopcomm.slave_comm, &local_n0, &local_0_start);


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

	// ==================================
	// Create and execute plan
	// ==================================
	plan = fftwl_mpi_plan_dft_r2c_2d(N0, N1, rho_x.data(), rho_k, loopcomm.slave_comm, FFTW_EXHAUSTIVE);
	fftwl_execute(plan);

	// ==================================
	// Divide by k^2
	// ==================================
	f0 = psi.dxdi * 2 * N0;
	f1 = psi.dydj * 2 * N1;
	for (long long int ind=0; i<rho_x_size; i++)
	{
		i = (ind + local_0_start) % N1;
		j = (ind + local_0_start) / N1;
		if (i <= N0/2)
		{
			kx = i / f0;
		} else {
			kx = (i-N0) / f0;
		}
		ky = j/f1;
		k2 = kx*kx + ky*ky;
		rho_k[ind][0] = rho_k[ind][0] / k2;
		rho_k[ind][1] = rho_k[ind][1] / k2;
	}

	// ==================================
	// Destroy plan
	// ==================================
	fftwl_destroy_plan(plan);

	// ==================================
	// Create and execute plan
	// ==================================
	real_out = fftwl_alloc_real(alloc_local);
	plan = fftwl_mpi_plan_dft_c2r_2d(N0, N1, rho_k, real_out, loopcomm.slave_comm, FFTW_EXHAUSTIVE);
	fftwl_execute(plan);

	// ==================================
	// Send Data
	// ==================================
	fftwl_send_local_size(local_n0, local_0_start);
	MPI_Send(real_out, rho_x_size, MPI_LONG_DOUBLE, 0, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);

	fftwl_free(rho_k);
	fftwl_destroy_plan(plan);

	// ==================================
	// Send wisdom back to master
	// ==================================
	fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);

	return 0;
}
