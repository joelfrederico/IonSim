#ifndef __SLAVE_LOOPS_H_INCLUDED__
#define __SLAVE_LOOPS_H_INCLUDED__

#include "simparams.h"
#include "scalar_data.h"
#include "ebeam.h"
#include "ions.h"

template<typename T>
int SL_get_rho(const unsigned int step_buf, const unsigned int substep_buf, const SimParams &simparams, ScalarData<T> &rho, const Ebeam &ebeam);

int SL_dump_ions(const SimParams &simparams, const LoopComm loopcomm, const unsigned int &step_buf, const unsigned int &substep_buf, const Ions &ions);

int SL_dump_electrons(const SimParams &simparams, const LoopComm loopcomm, const unsigned int &step_buf, const Ebeam &ebeam);

int SL_push_ions(const SimParams &simparams, const unsigned int &substep_buf, const Field_Data *field, Ions &ions);

int SL_get_efield(const SimParams &simparams, const Ebeam &ebeam, Field_Data *&field);

int SL_send_efield(const SimParams &simparams, Field_Data *&field);

int SL_get_ifield(const Field_Data *ion_field);

#endif
