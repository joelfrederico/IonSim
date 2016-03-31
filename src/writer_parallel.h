#ifndef __WRITER_PARALLEL_H_INCLUDED__
#define __WRITER_PARALLEL_H_INCLUDED__

#include "writer_base.h"
#include <mpi.h>

class WriterParallel : public WriterBase
{
	private:
		MPI::Intracomm *_comm_id_ptr;

	public:
		WriterParallel(const std::string &filename, MPI::Intracomm *comm_id);

};

#endif
