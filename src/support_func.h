#ifndef __SUPPORT_FUNC_H_INCLUDED__
#define __SUPPORT_FUNC_H_INCLUDED__

#include "mpi.h"
#include <string>
#include "classes.h"
#include "hdf5.h"

namespace ionsim
{
	double gamma2GeV(double gamma);
	double GeV2gamma(double GeV);
	double gaussian(double mean, double sigma);

	int sendloop(const int * message);
	int sendloop(const int * message, int step);

	int overwrite_file(std::string const &filename, MPI::Intracomm &slave_comm_id);
	hid_t open_file(std::string const &filename, MPI::Intracomm &slave_comm_id);
	int writeattribute(std::string const &attr_name, long attr_value, std::string const &filename, MPI::Intracomm &slave_comm_id);
	int writeattribute(std::string const &attr_name, int attr_value, std::string const &filename, MPI::Intracomm &slave_comm_id);
	int writeattribute(std::string const &attr_name, double attr_value, std::string const &filename, MPI::Intracomm &slave_comm_id);

	int dump(std::string const &filename, std::string const &group, std::string const &dataset, MPI::Intracomm &comm, Parts * ebeam);
	
	const int TAG_LOOP_INSTRUCT = 100;
	const int LOOP_KILL         = 1;
	const int LOOP_DUMP_IONS    = 2;
	const int LOOP_DUMP_E       = 3;
	const int LOOP_PUSH_E       = 4;
	const int LOOP_PUSH_IONS    = 5;
}

#endif
