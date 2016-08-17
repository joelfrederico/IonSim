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

	/* if ((j == 127) || (j == 128)) */
	/* { */
	/* 	return 1; */
	/* } else { */
	/* 	return 0; */
	/* } */

	/* long double i_mag, i_min, i_max, j_mag, j_min, j_max, r2; */

	/// return cos(valx/10)*cos(valy/10);

	/// return 0.5*cos(valy/10)+0.5;

	/// if ( (i % 4 == 0) || (j % 4 == 0))
	/// {
	/// 	return 1;
	/// } else {
	/// 	return 0;
	/// }

	/// if (valx*valx + valy*valy < 40)
	/// {
	/// 	return 1;
	/// } else {
	/// 	return 0;
	/// }

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
	N0 = 256;
	N1 = 256;
	long double *r_buf;
	fftwl_complex *c_buf;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	ptrdiff_t i_nonlocal;
	fftwl_plan planforward, planreverse;
	std::vector<std::complex<long double>> cdata;
	/* std::vector<long double> rdata; */
	long ind, ind2;
	long double delx, dely;
	delx = dely = 0.05;
	long double kx, ky, k2, fx, fy;
	long double val;

	alloc_local = fftwl_mpi_local_size_2d(N0, int(N1/2+1), loopcomm.slave_comm, &local_n0, &local_0_start);

	c_buf = fftwl_alloc_complex(alloc_local);
	r_buf = fftwl_alloc_real(2*alloc_local);

	planforward = fftwl_mpi_plan_dft_r2c_2d(N0, N1, r_buf, c_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE);

	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		i_nonlocal = i+local_0_start;
		for (ptrdiff_t j=0; j<N1; j++)
		{
			ind = i*(N1+2)+j;
			r_buf[ind] = f(i_nonlocal, j);
		}
	}

	fftwl_execute(planforward);

	cdata = MPI_convert_complex_vec(c_buf, local_n0, N1);

	MPI_Complex_div_k2(cdata, delx, dely, local_n0, local_0_start, N0, N1);

	MPI_Send_Scalar_Complex(cdata, local_n0, local_0_start, N0, N1);

	planreverse = fftwl_mpi_plan_dft_c2r_2d(N0, N1, c_buf, r_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE);

	MPI_convert_complex_buf(cdata, c_buf, local_n0, local_0_start, N1);

	fftwl_execute(planreverse);

	auto rdata = MPI_convert_real_vec(r_buf, local_n0, N0, N1);

	MPI_Send_Scalar_Real(rdata, local_n0, local_0_start, N0, N1);

	fftwl_free(r_buf);
	fftwl_free(c_buf);

	return 0;
}
