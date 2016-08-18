#ifndef __CONSTS_H_INCLUDED__
#define __CONSTS_H_INCLUDED__

#define JTF_PRINTVAL(x) std::cout << #x << ":\t" << x << std::endl
#define JTF_PRINTVAL_NOEND(x) std::cout << #x << ":\t" << x
#define JTF_PRINT(x) std::cout << #x << std::endl
#define JTF_PRINT_NOEND(x) std::cout << #x

#include <gsl/gsl_const_mksa.h>
#include <vector>
#include <complex>

// ========================================
// Define types
// ========================================
typedef long double ldouble;
typedef long long llong;
typedef std::vector<ldouble> ldouble_vec;
typedef int loopflag_t;
typedef int tag_t;
typedef int pushmethod_t;
typedef int parttype_t;
typedef std::complex<double> complex_double;
typedef int zdist_t;

// ========================================
// Define classes
// ========================================
const pushmethod_t PUSH_RUNGE_KUTTA = 1;
const pushmethod_t PUSH_SIMPLE      = 2;
const pushmethod_t PUSH_FIELD       = 3;

class PushMethod
{
	private:
		const int _PUSH_TYPE;

	public:
		PushMethod(const pushmethod_t PUSH_TYPE);
		bool operator==(const PushMethod &other) const;

		const std::string name;
};

// ========================================
// Loop flags
// ========================================
const zdist_t Z_DIST_FLAT  = 0;
const zdist_t Z_DIST_GAUSS = 1;

class zDist
{
	private:
		const zdist_t _Z_DIST;

	public:
		zDist(const zdist_t Z_DIST);

		const std::string name;
};

// ========================================
// Loop flags
// ========================================

const loopflag_t LOOP_KILL         = 1;
const loopflag_t LOOP_DUMP_IONS    = 2;
const loopflag_t LOOP_DUMP_E       = 3;
const loopflag_t LOOP_PUSH_E       = 4;
const loopflag_t LOOP_PUSH_IONS    = 5;
const loopflag_t LOOP_GET_EFIELD   = 6;
const loopflag_t LOOP_SEND_EFIELD  = 7;
const loopflag_t LOOP_GET_IFIELD   = 8;
const loopflag_t LOOP_RESET_IFIELD = 9;
const loopflag_t LOOP_FFTW_2D_DFT  = 10;
const loopflag_t LOOP_FFTW_2D_IDFT = 11;
const loopflag_t LOOP_GET_RHO      = 12;
const loopflag_t LOOP_START_E_ITER = 13;
const loopflag_t LOOP_START_I_ITER = 14;
const loopflag_t LOOP_GET_FIELDS   = 15;
const loopflag_t LOOP_GET_WISDOM   = 16;

// ========================================
// MPI tags
// ========================================
const tag_t TAG_LOOP_INSTRUCT       = 100;
const tag_t TAG_LOOP_MESSAGE        = 200;
const tag_t TAG_FIELD               = 300;
const tag_t TAG_COMPLEX_VEC_START   = 350;
const tag_t TAG_LDOUBLE_COMPLEX_VEC = 400;
const tag_t TAG_DOUBLE_COMPLEX_VEC  = 500;
const tag_t TAG_FLOAT_COMPLEX_VEC   = 600;
const tag_t TAG_LONG_LONG           = 700;
const tag_t TAG_MASTER_SLAVE        = 800;

const int MASTER_RANK = 0;

// ========================================
// Science stuff
// ========================================
const double ELECTRON_REST_ENERGY = GSL_CONST_MKSA_MASS_ELECTRON * GSL_CONST_MKSA_SPEED_OF_LIGHT * GSL_CONST_MKSA_SPEED_OF_LIGHT;

// ========================================
// Particle types
// ========================================
const parttype_t PARTS_E   = 2;
const parttype_t PARTS_ION = 1;

// ========================================
// HDF5 consts
// ========================================
// Number of entries to write at once
// because HDF5 crashes if you try to
// do too much all at once.
const int MAX_N_WRITE = 1e5;
#endif
