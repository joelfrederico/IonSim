#include "writer_parallel.h"
#include <mpi.h>

WriterParallel::WriterParallel(const std::string &filename, MPI::Intracomm *comm_id) : WriterBase(filename)
{
	_comm_id_ptr = comm_id;
}
