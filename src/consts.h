#ifndef __CONSTS_H_INCLUDED__
#define __CONSTS_H_INCLUDED__

#include <gsl/gsl_const_mksa.h>
#include <vector>
#include <complex>

typedef std::vector<double> double_vec;
typedef int loopflag_t;
typedef int tag_t;
typedef int pushmethod_t;
typedef int parttype_t;
typedef std::complex<double> complex_double;

const loopflag_t LOOP_KILL        = 1;
const loopflag_t LOOP_DUMP_IONS   = 2;
const loopflag_t LOOP_DUMP_E      = 3;
const loopflag_t LOOP_PUSH_E      = 4;
const loopflag_t LOOP_PUSH_IONS   = 5;
const loopflag_t LOOP_GET_EFIELD  = 6;
const loopflag_t LOOP_SEND_EFIELD = 7;

const tag_t TAG_LOOP_INSTRUCT = 100;
const tag_t TAG_LOOP_MESSAGE  = 200;
const tag_t TAG_FIELD         = 300;

const int MASTER_RANK = 0;

const double ELECTRON_REST_ENERGY = GSL_CONST_MKSA_MASS_ELECTRON * GSL_CONST_MKSA_SPEED_OF_LIGHT * GSL_CONST_MKSA_SPEED_OF_LIGHT;

const parttype_t PARTS_E   = 2;
const parttype_t PARTS_ION = 1;

const pushmethod_t PUSH_RUNGE_KUTTA = 1;
const pushmethod_t PUSH_SIMPLE      = 2;
const pushmethod_t PUSH_FIELD       = 3;

// Number of entries to write at once
// because HDF5 crashes if you try to
// do too much all at once.
const int MAX_N_WRITE = 1e5;
#endif
