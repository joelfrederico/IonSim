#ifndef __SUPPORT_FUNC_H_INCLUDED__
#define __SUPPORT_FUNC_H_INCLUDED__

#include <mpi.h>
#include <string>
#include "parts.h"
#include <hdf5.h>
#include <gsl/gsl_const_mksa.h>
#include "field_data.h"
#include "field_comm.h"

namespace ionsim
{
	// ==================================
	// Methods
	// ==================================
	double GeV2gamma(double GeV);
	double gamma2GeV(double gamma);
	double gaussian(double mean, double sigma);
	int sendloop(const int &message);
	int loop_get_fields(Field_Comm &fieldcomm, Field_Data &field);
	/* int loop_push_ions(Field &field); */
	int sendloop(const int &message, int step);

	// ==================================
	// Consts
	// ==================================
	const double ELECTRON_REST_ENERGY = GSL_CONST_MKSA_MASS_ELECTRON * pow(GSL_CONST_MKSA_SPEED_OF_LIGHT, 2);

	const int TAG_LOOP_INSTRUCT = 100;
	const int TAG_FIELD         = 200;

	const loopflag_t LOOP_DUMP_E     = 3;
	const loopflag_t LOOP_DUMP_IONS  = 2;
	const loopflag_t LOOP_GET_EFIELD = 6;
	const loopflag_t LOOP_KILL       = 1;
	const loopflag_t LOOP_PUSH_E     = 4;
	const loopflag_t LOOP_PUSH_IONS  = 5;

	const parttype_t PARTS_E   = 2;
	const parttype_t PARTS_ION = 1;

	const pushmethod_t PUSH_RUNGE_KUTTA = 1;
	const pushmethod_t PUSH_SIMPLE      = 2;
	const pushmethod_t PUSH_FIELD       = 3;

	const int MAX_N_WRITE = 1e5;
}

#endif
