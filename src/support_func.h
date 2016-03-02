#ifndef __SUPPORT_FUNC_H_INCLUDED__
#define __SUPPORT_FUNC_H_INCLUDED__

#include "mpi.h"
#include <string>
#include "classes.h"

namespace ionsim
{
	double gamma2GeV(double gamma);
	double GeV2gamma(double GeV);
	int dump(std::string const &filename, MPI::Intracomm &comm, Ebeam *ebeam);
}

#endif
