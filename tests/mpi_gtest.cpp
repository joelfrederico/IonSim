#include <mpi.h>
#include <hdf5.h>

int main(int argc, char **argv)
{
	MPI_Info info = MPI_INFO_NULL;
	MPI_Status status;
	MPI_Comm comm = MPI_COMM_WORLD;
	int size;
	int rank;
	hid_t fapl, fh;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	std::cout << "Rank: " << rank << std::endl;

	fapl = H5Pcreate(H5P_FILE_ACCESS);

	H5Pset_fapl_mpio(fapl, comm, info);

	fh = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

	H5Pclose(fapl);
	H5Fclose(fh);

	MPI_Finalize();

	return 0;
}
