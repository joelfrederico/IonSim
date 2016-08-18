#include "mpi_vec.h"
#include "loop_comm.h"
#include "scalar_data.h"
#include <fftw3-mpi.h>
#include "writer_serial.h"
#include "support_func.h"

// ========================================
// Master: Master(scalar<real>)->Slave(cbuf)
// ========================================
template<typename T>
int MPI_Master_Send_Scalar_real_to_buf(const ScalarData<T> &buf)
{
	LoopComm loopcomm;
	MPI_Datatype mpi_type;
	ptrdiff_t local_n0, local_0_start, N0, N1;
	T *tbuf;

	JTF_PRINT_NOEND(Rank: ) << buf.x_pts_vec().size() << std::endl;;

	mpi_type = ionsim::convert_typeid_to_mpi<T>();

	for (int id=1; id<loopcomm.p; id++)
	{
		MPI_Recv_local(id, local_n0, local_0_start, N0, N1);
		tbuf = buf.vdata().data() + local_0_start;
		MPI_Send(tbuf, local_n0*N1, mpi_type, id, TAG_MASTER_SLAVE, MPI_COMM_WORLD);
	}
	return 0;
}
template int MPI_Master_Send_Scalar_real_to_buf(const ScalarData<long double> &buf);

// ========================================
// Slave: Master(scalar<real>)->Slave(cbuf)
// ========================================
template<typename T>
int MPI_Slave_Recv_buf_from_Scalar_real(T *buf, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1)
{
	MPI_Status status;
	auto mpitype = ionsim::convert_typeid_to_mpi<T>();
	int count;

	MPI_Send_local(local_n0, local_0_start, N0, N1);

	MPI_Probe(0, TAG_MASTER_SLAVE, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, mpitype, &count);

	if (count != local_n0*N1) throw std::runtime_error("Receiving incommensurate number of counts!");

	MPI_Recv(buf, local_n0*N1, mpitype, 0, TAG_MASTER_SLAVE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	return 0;
}
template int MPI_Slave_Recv_buf_from_Scalar_real(long double *buf, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

int MPI_Send_local(const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1)
{
	std::vector<long long> buf = {local_n0, local_0_start, N0, N1};
	MPI_Send(buf.data(), buf.size(), MPI_LONG_LONG_INT, 0, TAG_LONG_LONG, MPI_COMM_WORLD);
	return 0;
}

int MPI_Recv_local(const int id, ptrdiff_t &local_n0, ptrdiff_t &local_0_start, ptrdiff_t &N0, ptrdiff_t &N1)
{
	std::vector<long long> buf(4);

	MPI_Recv(buf.data(), buf.size(), MPI_LONG_LONG_INT, id, TAG_LONG_LONG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	local_n0      = buf[0];
	local_0_start = buf[1];
	N0            = buf[2];
	N1            = buf[3];
	return 0;
}

template<typename T>
int MPI_Send_vector(const std::vector<T> &buf)
{
	MPI_Datatype mpitype;
	int tag;
	int count = buf.size();

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

	if (buf.size() < count) throw std::runtime_error("Vector smaller than count.");
	MPI_Send(buf.data(), count, mpitype, 0, tag, MPI_COMM_WORLD);

	return 0;
}

template<typename T>
int MPI_Recv_vector(int id, std::vector<T> &buf)
{
	MPI_Status status;
	MPI_Datatype mpitype;
	int count;
	long ind;

	MPI_Probe(id, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	if (status.MPI_TAG == TAG_LDOUBLE_COMPLEX_VEC)
	{
		if (typeid(T) != typeid(long double)) throw std::runtime_error("Types don't match");
		mpitype = MPI_LONG_DOUBLE;
	} else if (status.MPI_TAG == TAG_DOUBLE_COMPLEX_VEC) {
		if (typeid(T) != typeid(double)) throw std::runtime_error("Types don't match");
		mpitype = MPI_DOUBLE;
	} else if (status.MPI_TAG == TAG_FLOAT_COMPLEX_VEC) {
		if (typeid(T) != typeid(float)) throw std::runtime_error("Types don't match");
		mpitype = MPI_FLOAT;
	}

	MPI_Get_count(&status, mpitype, &count);
	buf.resize(count);

	if (buf.size() < count) throw std::runtime_error("Vector smaller than count.");
	MPI_Recv(buf.data(), count, mpitype, id, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	return 0;
}

template<typename T>
int MPI_Send_complex(const std::vector<std::complex<T>> &buf)
{
	LoopComm loopcomm;
	typename std::vector<T> dbuf;
	long ind;
	int ccount;

	ccount = buf.size();
	dbuf.resize(2*ccount);

	for (long i=0; i<ccount; i++)
	{
		ind = i*2;
		if (i > ccount) throw std::runtime_error("Index greater than size.");
		if (ind+1 > 2*ccount) throw std::runtime_error("Index greater than size.");
		dbuf[ind]   = buf[i].real();
		dbuf[ind+1] = buf[i].imag();
	}

	/* std::cout << "Sending size: " << ccount << std::endl; */
	return MPI_Send_vector(dbuf);
}

template<typename T>
int MPI_Recv_complex(int id, std::vector<std::complex<T>> &buf)
{
	int ccount;
	long ind;
	typename std::vector<T> dbuf;

	MPI_Recv_vector(id, dbuf);

	ccount = dbuf.size()/2;
	/* std::cout << "Received size: " << ccount << std::endl; */

	buf.resize(ccount);
	for (long i=0; i<ccount; i++)
	{
		ind = i*2;
		buf[i].real(dbuf[ind]);
		buf[i].imag(dbuf[ind+1]);
	}

	return 0;
}

template<typename T>
int MPI_Recv_Scalar_Real(ScalarData<T> &rdata)
{
	LoopComm loopcomm;
	std::vector<long double> rbuf;
	typename decltype(rdata.vdata())::size_type i_nonlocal;
	ptrdiff_t local_n0, local_0_start, N0, N1;

	for (int id=1; id < loopcomm.p; id++)
	{
		MPI_Recv_local(id, local_n0, local_0_start, N0, N1);
		MPI_Recv_vector(id, rbuf);
		for (long i=0; i<local_n0; i++)
		{
			i_nonlocal = i+local_0_start;
			for (decltype(i_nonlocal) j=0; j<N1; j++)
			{
				rdata.ind(i_nonlocal, j) = rbuf[i*N1+j];
			}
		}
	}
	return 0;
}

template<typename T>
int MPI_Recv_Scalar_Complex(ScalarData<std::complex<T>> &cdata)
{
	LoopComm loopcomm;
	typename std::vector<std::complex<T>> cbuf;
	typename decltype(cdata.vdata())::size_type i_nonlocal;
	ptrdiff_t local_n0, local_0_start, N0, N1;
	for (int id=1; id < loopcomm.p; id++)
	{
		MPI_Recv_local(id, local_n0, local_0_start, N0, N1);
		if ((cdata.x_pts(0) != N0) || cdata.x_pts(1) != (N1/2+1)) throw std::runtime_error("Sizes not commensurate");
		MPI_Recv_complex(id, cbuf);
		for (long i=0; i<local_n0; i++)
		{
			i_nonlocal = i + local_0_start;
			for (decltype(i_nonlocal) j=0; j<(N1/2+1); j++)
			{
				cdata.ind(i_nonlocal, j) = cbuf[i*(N1/2+1)+j];
			}
		}
	}

	auto j = i_nonlocal;
	i_nonlocal = 1;
	j=0;

	JTF_PRINTVAL(cdata.ind(i_nonlocal, j).real());
	JTF_PRINTVAL(cdata.ind(i_nonlocal, j).imag());
	return 0;
}

template<typename T>
int MPI_Send_Scalar_Real(const std::vector<T> &rdata, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1)
{
	MPI_Send_local(local_n0, local_0_start, N0, N1);
	MPI_Send_vector(rdata);
	return 0;
}

std::vector<long double> MPI_convert_real_vec(long double *r_buf, const ptrdiff_t local_n0, const ptrdiff_t N0, const ptrdiff_t N1)
{
	long ind, ind2, maxallowed;
	maxallowed = local_n0*N1;
	std::vector<long double> rdata;

	rdata.resize(maxallowed);
	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		for (ptrdiff_t j=0; j<N1; j++)
		{
			ind = i*(2*(N1/2+1))+j;
			ind2 = i*N1+j;
			rdata[ind2] = r_buf[ind] / (N0*N1);
		}
	}
	return rdata;
}

std::vector<std::complex<long double>> MPI_convert_complex_vec(fftwl_complex *c_buf, ptrdiff_t local_n0, ptrdiff_t N1)
{
	typename std::vector<std::complex<long double>> cdata;
	long csize;
	
	csize = local_n0*(N1/2+1);
	cdata.resize(csize);
	for (long i=0; i<csize; i++)
	{
		cdata[i].real(c_buf[i][0]);
		cdata[i].imag(c_buf[i][1]);
	}
	return cdata;
}

template<typename T>
int MPI_Send_Scalar_Complex(const std::vector<std::complex<T>> cdata, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1)
{
	MPI_Send_local(local_n0, local_0_start, N0, N1);
	MPI_Send_complex(cdata);
	return 0;
}

template<typename T>
int MPI_Complex_div_k2(std::vector<std::complex<T>> &cdata, const T delx, const T dely, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1)
{
	T fx, fy, kx, ky, k2;
	ptrdiff_t i_nonlocal;
	long ind;

	fx = delx*2;
	fy = dely*2;
	for (long i=0; i<local_n0; i++)
	{
		i_nonlocal = i+local_0_start;
		if (i_nonlocal <= N0/2)
		{
			kx = i_nonlocal;
		} else {
			kx = i_nonlocal-N0;
		}
		kx /= fx;

		for (long j=0; j<(N1/2+1); j++)
		{
			ind = i*(N1/2+1) + j;

			ky = j;
			ky /= fy;
			
			k2 = kx*kx + ky*ky;
			if ((i_nonlocal == 0) && (j==0))
			{
				cdata[ind] = 0;
			} else {
				cdata[ind] /= k2;
			}

			/* if ((i_nonlocal == 1) && (j == 0)) */
			/* { */
			/* 	JTF_PRINTVAL(local_0_start); */
			/* 	JTF_PRINTVAL(i_nonlocal); */
			/* 	JTF_PRINTVAL(j); */
			/* 	JTF_PRINTVAL(ind); */
			/* 	JTF_PRINTVAL(N0); */
			/* 	JTF_PRINTVAL(N1); */
			/* 	JTF_PRINTVAL(N1/2+1); */
			/* 	cdata[ind].real(100); */
			/* 	cdata[ind].imag(200); */
			/* } */
		}
	}

	return 0;
}

int MPI_convert_complex_buf(const std::vector<std::complex<long double>> cdata, fftwl_complex *c_buf, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N1)
{
	ptrdiff_t i_nonlocal;
	long ind;
	long maxallowed = cdata.size();

	for (ptrdiff_t i=0; i<local_n0; i++)
	{
		i_nonlocal=i+local_0_start;
		for(ptrdiff_t j=0; j<(N1/2+1); j++)
		{
			ind = i*(N1/2+1)+j;
			if (ind >= maxallowed) JTF_PRINT(Over max);
			c_buf[ind][0] = cdata[ind].real();
			c_buf[ind][1] = cdata[ind].imag();
		}
	}
	return 0;
}


// ========================================
// Instantiation
// ========================================

template int MPI_Send_complex(const std::vector<std::complex<long double>> &buf);

template int MPI_Recv_complex(int id, std::vector<std::complex<long double>> &buf);

template int MPI_Recv_Scalar_Real(ScalarData<long double> &rdata);

template int MPI_Send_Scalar_Real(const std::vector<long double> &rdata, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

template int MPI_Send_Scalar_Complex(std::vector<std::complex<long double>> cdata, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);

template int MPI_Recv_Scalar_Complex(ScalarData<std::complex<long double>> &cdata);

template int MPI_Complex_div_k2(std::vector<std::complex<long double>> &cdata, const long double delx, const long double dely, const ptrdiff_t local_n0, const ptrdiff_t local_0_start, const ptrdiff_t N0, const ptrdiff_t N1);
