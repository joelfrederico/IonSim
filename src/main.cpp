#include "config.h"
#include "main.h"

#include "master.h"
#include "slave.h"

#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <gflags/gflags.h> // #include <google/gflags.h>
#include <fftw3-mpi.h>

const int INT_TAG = 1;

DEFINE_bool(verbose, true, "Verbose");
DEFINE_string(file, "config.xml", "Configuration file");
DEFINE_string(wisdom_file, ".fftw-wisdom", "FFTW wisdom file");

int main(int argc, char **argv)
{
	gflags::ParseCommandLineFlags(&argc, &argv, true);

	int id, c;
	bool verbose = FLAGS_verbose;
	
	// ==============================
	// Initialize MPI
	// ==============================
	MPI_Init(&argc, &argv);
	fftwl_mpi_init();

	// ==============================
	// Get the individual process ID.
	// ==============================
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if (verbose && id==0) printf("OpenMPI v%d.%d.%d\n", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION);

	// ==============================
	// Simulate things
	// ==============================
	if (id == 0) {
		master();
	} else {
		slave();
	}

	MPI_Barrier(MPI_COMM_WORLD);

	fftwl_mpi_cleanup();
	MPI_Finalize();

	if (verbose)
	{
		if (id == 0)
		{
			printf("Master finished\n");
		}
		else
		{
			printf("Slave %d finished\n", id);
		}
	}

	return 0;
}
