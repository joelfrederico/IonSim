#ifndef __MASTER_LOOP_H_INCLUDED__
#define __MASTER_LOOP_H_INCLUDED__

#include "field_data.h"
#include "writer_serial.h"
#include "simparams.h"

int ML_overwrite_file(const SimParams &simparams);

int ML_write_field(const SimParams &simparams, const unsigned int e_step, const Field_Data &field);

template<typename T>
int ML_write_scalar(const SimParams &simparams, const ScalarData<T> &rho)
{
	WriterSerial writer_s(simparams.filename);
	writer_s.writedata(rho, "rho");
}

#endif
