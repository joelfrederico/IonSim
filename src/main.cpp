#include "config.h"
#include "main.h"

#include "master.h"
#include "slave.h"

#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <gflags/gflags.h> // #include <google/gflags.h>

const int INT_TAG = 1;

DEFINE_bool(verbose, true, "Verbose");
DEFINE_string(file, "config.xml", "Configuration file");

int main(int argc, char **argv)
{
	gflags::ParseCommandLineFlags(&argc, &argv, true);

	int id, c;
	bool verbose = FLAGS_verbose;
	
	// ==============================
	// Initialize MPI
	// ==============================
	MPI_Init(&argc, &argv);

	// ==============================
	// Get the individual process ID.
	// ==============================
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if (verbose && id==0) printf("OpenMPI v%d.%d.%d\n", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION);

	// ==============================
	// Simulate things
	// ==============================
	if (id == 0) {
		master(verbose);
	} else {
		slave(verbose);
	}

	MPI_Barrier(MPI_COMM_WORLD);

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
