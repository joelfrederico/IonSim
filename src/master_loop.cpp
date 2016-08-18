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

template<typename T>
int ML_SolvePoisson(const SimParams &simparams, const std::string wisdom_file, const ScalarData<T> &rho, ScalarData<T> &psi)
{
	decltype(rho.x_pts_vec()) x_pts_complex = rho.x_pts_vec();
	x_pts_complex[1] = x_pts_complex[1]/2 + 1;
	auto edge_mag_complex = rho.edge_mag_vec();

	ScalarData<std::complex<long double>> cdata(x_pts_complex, edge_mag_complex);
	auto rdata = rho;

	MPI_Recv_Scalar_Complex(cdata);
	MPI_Recv_Scalar_Real(rdata);

	{
		WriterSerial writer_s(simparams.filename);
		writer_s.writedata(cdata, "complex");
	}

	{
		WriterSerial writer_s(simparams.filename);
		writer_s.writedata(rdata, "rdata");
	}

	fftwl_mpi_gather_wisdom(MPI_COMM_WORLD);
	fftwl_export_wisdom_to_filename(wisdom_file.c_str());

	return 0;
}

template int ML_SolvePoisson(const SimParams &simparams, const std::string wisdom_file, const ScalarData<long double> &rho, ScalarData<long double> &psi);
