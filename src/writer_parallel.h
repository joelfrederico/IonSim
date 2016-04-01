#ifndef __WRITER_PARALLEL_H_INCLUDED__
#define __WRITER_PARALLEL_H_INCLUDED__

#include "writer_base.h"
#include <mpi.h>
#include "parts.h"

class WriterParallel : public WriterBase
{
	private:
		int p;
		int id;
		MPI::Intracomm *comm_id_ptr;
		int _init(const std::string &filename, MPI::Intracomm *comm_id, bool overwrite);

	public:
		WriterParallel(const std::string &filename, MPI::Intracomm *comm_id);
		WriterParallel(const std::string &filename, MPI::Intracomm *comm_id, bool overwrite);

		int write_attributes(const SimParams &simparams) const;
		int open_file(std::string const &filename);
		int overwrite_file_parallel(const std::string &filename);
		int writedata(long step, std::string const &dataset, const Parts &parts);

};

#endif
