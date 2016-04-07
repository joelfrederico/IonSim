#ifndef __WRITER_PARALLEL_H_INCLUDED__
#define __WRITER_PARALLEL_H_INCLUDED__

#include "writer_base.h"
#include <mpi.h>
#include "parts.h"
#include "loop_comm.h"

class WriterParallel : public WriterBase
{
	private:
		const LoopComm loopcomm;
		const MPI_Comm _comm_id;

		int _init(const std::string &filename, bool overwrite);
		int _writedata(hid_t &group_id, const std::string &dataset, const Parts &parts);

	public:
		WriterParallel(const std::string &filename, const MPI_Comm comm_id);
		WriterParallel(const std::string &filename, const MPI_Comm comm_id, bool overwrite);

		int write_attributes(const SimParams &simparams) const;
		int open_file();
		int overwrite_file_parallel();
		int writedata(long step, const std::string &dataset, const Parts &parts);
		int writedata_substep(long step, long substep, const std::string &dataset, const std::string &subgroup, const Parts &parts);

};

#endif
