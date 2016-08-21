#ifndef __MPI_VEC__H_INCLUDED__
#define __MPI_VEC__H_INCLUDED__

#include <vector>
#include <complex>
#include "scalar_data.h"
#include <fftw3-mpi.h>

template<typename T>
int MPI_Master_Send_Scalar_real_to_buf(ScalarData<T> &buf);

template<typename T>
int MPI_Slave_Recv_buf_from_Scalar_real(T *buf, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

int MPI_Send_local(const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

int MPI_Recv_local(const int id, ptrdiff_t &local_n0, ptrdiff_t &local_0_start, ptrdiff_t &N0, ptrdiff_t &N1);

template<typename T>
int MPI_Send_complex(const std::vector<std::complex<T>> &buf);

template<typename T>
int MPI_Recv_complex(int id, std::vector<std::complex<T>> &buf);

template<typename T>
int MPI_Send_vector(const std::vector<T> &buf);

template<typename T>
int MPI_Recv_vector(int id, std::vector<T> &buf);

template<typename T>
int MPI_Recv_Scalar_Real(ScalarData<T> &rdata);

template<typename T>
int MPI_Recv_Scalar_Complex(ScalarData<std::complex<T>> &cdata);

template<typename T>
int MPI_Send_Scalar_Real(const std::vector<T> &rdata, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

template<typename T>
int MPI_Send_Scalar_Complex(const std::vector<std::complex<T>> cdata, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

template<typename T>
int MPI_Complex_div_k2(std::vector<std::complex<T>> &cdata, const T delx, const T dely, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

int SolvePoisson();

// ==================================
// Don't Need Instantiation
// ==================================
std::vector<long double> MPI_convert_real_vec(long double *r_buf, const ptrdiff_t local_n0, const ptrdiff_t N0, const ptrdiff_t N1);
std::vector<std::complex<long double>> MPI_convert_complex_vec(fftwl_complex *c_buf, const ptrdiff_t local_n0, const ptrdiff_t N1);
int MPI_convert_complex_buf(const std::vector<std::complex<long double>> cdata, fftwl_complex *c_buf, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N1);


#endif
