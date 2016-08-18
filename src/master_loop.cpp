#include "master_loop.h"
#include "writer_serial.h"
#include "consts.h"
#include "fftw_classes.h"
#include "mpi_vec.h"

int ML_overwrite_file(const SimParams &simparams)
{
	WriterSerial writer_s(simparams.filename, true);
	writer_s.write_attributes(simparams);

	return 0;
}

int ML_write_field(const SimParams &simparams, const unsigned int e_step, const Field_Data &field)
{
	WriterSerial writer_s(simparams.filename);
	writer_s.writedata(e_step, field);

	return 0;
}

template<typename T>
int ML_loop_get_efields(const LoopComm &loopcomm, const unsigned int e_step, const ScalarData<T> &rho)
{
	long long local_n0, local_0_start;
	long long rho_size;
	std::vector<long double> buf;
	long long i_nonlocal;
	ScalarData<std::complex<T>> cdata(rho.x_pts, rho.y_pts/2+1, 1, 1, 1, 1);

	loopcomm.instruct(LOOP_GET_FIELDS);
	for (int id=1; id < loopcomm.p; id++)
	{
		fftwl_recv_local_size(local_n0, local_0_start, id);

		rho_size = local_n0*rho.y_pts;

		buf.resize(rho_size);
		for (unsigned long i=0; i<local_n0; i++)
		{
			i_nonlocal = i + local_0_start;
			for (unsigned long j=0; j<rho.y_pts; j++)
			{
				buf[ionsim::row_major(i, j, rho.y_pts)] = rho.ind(i_nonlocal, j, e_step);
			}
		}

		// Send rho
		MPI_Send(buf.data(), rho_size, MPI_LONG_DOUBLE, id, TAG_LOOP_MESSAGE, MPI_COMM_WORLD);
	}

	for (int id=1; id < loopcomm.p; id++)
	{
		MPI_Recv_complex(id, cdata.vdata());
	}

	return 0;
}

int ML_SolvePoisson(const char *wisdom_file)
{
	std::vector<unsigned long> x_pts_complex = {256, 129};
	std::vector<unsigned long> x_pts_real    = {256, 256};
	std::vector<long double> edge_mag        = {1, 1};

	ScalarData<std::complex<long double>> cdata(x_pts_complex, edge_mag);
	ScalarData<long double> rdata(x_pts_real, edge_mag);
	std::vector<std::complex<long double>> cbuf;
	std::vector<long double> rbuf;
	WriterSerial *writer_s;
	ptrdiff_t local_n0, local_0_start, N0, N1;


	MPI_Recv_Scalar_Complex(cdata);
	MPI_Recv_Scalar_Real(rdata);

	writer_s = new WriterSerial("myfile.h5", true);
	writer_s->writedata(cdata, "complex");
	delete writer_s;

	writer_s = new WriterSerial("myfile.h5");
	writer_s->writedata(rdata, "rdata");
	delete writer_s;

	fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);
	fftwl_export_wisdom_to_filename(wisdom_file);

	return 0;
}
