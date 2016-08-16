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

long c_ind(const long i, const long j, const ptrdiff_t N1)
{
	return i*N1 + j;
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
		valx = 1;
	} else {
		valx = 0;
	}
	return valx;
}

int psifftw_base(LoopComm loopcomm)
{
	ptrdiff_t N0 = 256;
	ptrdiff_t N1 = 256;
	long double *r_buf;
	fftwl_complex *c_buf;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	ptrdiff_t i_nonlocal;
	fftwl_plan planforward, planreverse;
	std::vector<std::complex<long double>> cdata;
	std::vector<long double> rdata;
	long ind, ind2, maxallowed, csize;

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

	MPI_Send_Scalar_Complex(c_buf, local_n0, local_0_start, N0, N1);

	csize = local_n0*(N1/2+1);
	cdata.resize(csize);
	for (long i=0; i<csize; i++)
	{
		cdata[i].real(c_buf[i][0]);
		cdata[i].imag(c_buf[i][1]);
	}


	planreverse = fftwl_mpi_plan_dft_c2r_2d(N0, N1, c_buf, r_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE);

	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		i_nonlocal=i+local_0_start;
		for(ptrdiff_t j=0; j<(N1/2+1); j++)
		{
			ind = i*(N1/2+1)+j;
			if (ind > maxallowed) JTF_PRINT(Over max);
			c_buf[ind][0] = cdata[ind].real();
			c_buf[ind][1] = cdata[ind].imag();
		}
	}

	fftwl_execute(planreverse);

	maxallowed = local_n0*N1;
	rdata.resize(maxallowed);
	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		i_nonlocal = i+local_0_start;
		for (ptrdiff_t j=0; j<N1; j++)
		{
			ind = i*(2*(N1/2+1))+j;
			ind2 = i*N1+j;
			rdata[ind2] = r_buf[ind] / (N0*N1);
		}
	}

	MPI_Send_local(local_n0, local_0_start, N0, N1);
	MPI_Send_vector(rdata);

	fftwl_free(r_buf);
	fftwl_free(c_buf);

	return 0;
}

int lyingness(LoopComm loopcomm)
{
	// ==================================
	// Allocate variables
	// ==================================
	fftwl_plan planforward;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	/// ldouble_vec rho_x;
	long double val;
	/// long long int rho_x_size;
	/// long long int i, j;
	ptrdiff_t i_nonlocal;

	/// ScalarData<ldouble> psi(128, 128, 1, 1, 1, 1);
	/// long double f0, f1, kx, ky, k2;
	long double *r_buf;
	/// std::vector<long double> real_out;
	fftwl_complex *c_buf;
	/// long ind, ind2;
	std::vector<std::complex<long double>> cdata;

	// ==================================
	// Size of FFT
	// ==================================
	/* ptrdiff_t N0 = simparams.n_field_x; */
	/* ptrdiff_t N1 = simparams.n_field_y; */
	ptrdiff_t N0 = 128;
	ptrdiff_t N1 = 128;

	// ==================================
	// Determine local size
	// ==================================
	alloc_local = fftwl_mpi_local_size_2d(N0, N1/2+1, loopcomm.slave_comm, &local_n0, &local_0_start);

	// ==================================
	// Send local data
	// ==================================
	// fftwl_send_local_size(local_n0, local_0_start);
	// rho_x_size = local_n0*N1;
	// rho_x.resize(rho_x_size);
	// MPI_Recv(rho_x.data(), rho_x_size, MPI_LONG_DOUBLE, 0, TAG_LOOP_MESSAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// ==================================
	// Allocate complex data
	// LARGER THAN NECESSARY
	// ==================================
	/// rho_k = fftwl_alloc_complex(alloc_local);
	c_buf = fftwl_alloc_complex(alloc_local);
	r_buf = fftwl_alloc_real(2*alloc_local);

	/// real_out.resize(local_n0*N1);

	// JTF_PRINTVAL_NOEND(alloc_local)   << ", id: "   << loopcomm.id           << std::endl;
	// JTF_PRINTVAL_NOEND(2*alloc_local) << ", id: "   << loopcomm.id           << std::endl;
	// JTF_PRINTVAL_NOEND(local_n0)      << ", id: "   << loopcomm.id           << std::endl;
	// JTF_PRINTVAL_NOEND(local_0_start) << ", id: "   << loopcomm.id           << std::endl;
	// JTF_PRINTVAL_NOEND(r_buf)         << ", size: " << sizeof(long double)   << ", id: "   << loopcomm.id << std::endl;
	// JTF_PRINTVAL_NOEND(c_buf)         << ", size: " << sizeof(fftwl_complex) << ", id: "   << loopcomm.id << std::endl;
	/* JTF_PRINTVAL_NOEND(real_out.data()) << ", size: " << sizeof(long double)   << ", id: "   << loopcomm.id << std::endl; */
	/* JTF_PRINTVAL_NOEND(rho_k)           << ", size: " << sizeof(fftwl_complex) << ", id: "   << loopcomm.id << std::endl; */

	// ==================================
	// Create and execute plan
	// ==================================
	/// planforward = fftwl_mpi_plan_dft_r2c_2d(N0, N1, r_buf, c_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);

	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		i_nonlocal = i+local_0_start;
		for (ptrdiff_t j=0; j<N1; j++)
		{
			if ((i_nonlocal==64) && (j==64))
			{
				JTF_PRINTVAL(i);
				JTF_PRINTVAL(i_nonlocal);
				JTF_PRINTVAL(j);
				val = 1;
			} else {
				val = 0;
			}
			r_buf[i*(2*(N1/2+1))+j] = val;
		}
	}

	/* long max = 0; */
	/* for (long i=0; i < local_n0; i++) */
	/* { */
	/* 	for (long j=0; j < N1; j++) */
	/* 	{ */
	/* 		val = rho_x[ionsim::row_major(i, j, N1)]; */
	/* 		ind = i*(2*(N1/2+1))+j; */

	/* 		r_buf[ind] = val; */

	/* 		if (ind > max) max = ind; */
	/* 		if (val != 0) JTF_PRINTVAL(val); */
	/* 	} */
	/* } */

	/// fftwl_execute(planforward);

	cdata.resize(local_n0*(N1/2+1));
	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		i_nonlocal=i+local_0_start;
		for(ptrdiff_t j=0; j<N1; j++)
		{
			cdata[i*N1+j].real(c_buf[i*N1+j][0]);
			cdata[i*N1+j].imag(c_buf[i*N1+j][1]);
		}
	}

	// // ==================================
	// // Divide by k^2
	// // ==================================
	// f0 = psi.dxdi * 2 * N0;
	// f1 = psi.dydj * 2 * N1;
	// for (long i=0; i < local_n0; i++)
	// {
	// 	i_nonlocal = i + local_0_start;
	// 	if (i_nonlocal <= N0/2)
	// 	{
	// 		kx = i_nonlocal / f0;
	// 	} else {
	// 		kx = (i_nonlocal-N0) / f0;
	// 	}
	// 	for (long j=0; j< N1/2+1; j++)
	// 	{
	// 		ky = j/f1;
	// 		k2 = kx*kx + ky*ky;
	// 		ind = c_ind(i, j, N1/2+1);
	// 		if ((i_nonlocal==0) && (j==0))
	// 		{
	// 			rho_k[ind][0] = 0;
	// 			rho_k[ind][1] = 0;
	// 		} else {
	// 			rho_k[ind][0] = c_buf[ind][0] / k2;
	// 			rho_k[ind][1] = c_buf[ind][1] / k2;
	// 		}
	// 		ind2 = ionsim::row_major(i, j, N0);
	// 		cdata[ind2].real(rho_k[ind][0]);
	// 		cdata[ind2].imag(rho_k[ind][1]);
	// 	}
	// }

	/* cdata[0].real(10); */
	/* cdata[0].imag(-10); */
	/// MPI_Send_complex(cdata, local_n0, local_0_start);

	fftwl_free(r_buf);
	/// fftwl_free(rho_k);
	fftwl_free(c_buf);
	/// fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);
	return 0;

	// fftwl_free(c_buf);
	// fftwl_free(r_buf);

	// // ==================================
	// // Create and execute plan
	// // ==================================
	// alloc_local = fftwl_mpi_local_size_2d(N0, N1/2+1, loopcomm.slave_comm, &local_n0, &local_0_start);
	// c_buf = fftwl_alloc_complex(alloc_local);
	// r_buf = fftwl_alloc_real(2*alloc_local);

	// planreverse = fftwl_mpi_plan_dft_c2r_2d(N0, N1, c_buf, r_buf, loopcomm.slave_comm, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
	// JTF_PRINT(Or try dying);

	// for (long i=0; i < N0; i++)
	// {
	// 	for (long j=0; j<N1/2+1; j++)
	// 	{
	// 		c_buf[c_ind(i, j, N1/2+1)][0] = rho_k[ionsim::row_major(i, j, N1)][0];
	// 		c_buf[c_ind(i, j, N1/2+1)][1] = rho_k[ionsim::row_major(i, j, N1)][1];
	// 	}
	// }

	// JTF_PRINT(Yaya);

	// fftwl_execute(planreverse);

	// JTF_PRINT(Nana);
	// for (long i=0; i < local_n0; i++)
	// {
	// 	for (long j=0; j<N1; j++)
	// 	{
	// 		real_out[ionsim::col_major(i, j, local_n0)] = r_buf[ionsim::row_major(i, j, N1)];
	// 	}
	// }

	// JTF_PRINT(Further);

	// // ==================================
	// // Send Data
	// // ==================================
	// fftwl_send_local_size(local_n0, local_0_start);
	// MPI_Send(real_out.data(), rho_x_size, MPI_LONG_DOUBLE, 0, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);

	// fftwl_free(r_buf);
	// fftwl_free(rho_k);

	// // ==================================
	// // Send wisdom back to master
	// // ==================================
	// fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);

	// fftwl_destroy_plan(planforward);
	// fftwl_destroy_plan(planreverse);

	return 0;
}
