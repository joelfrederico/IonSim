#ifndef __SLAVE_H_INCLUDED__
#define __SLAVE_H_INCLUDED__

#include "mpi.h"

int slave(int &p, int &id, MPI::Intracomm &slave_comm_id);


#endif
