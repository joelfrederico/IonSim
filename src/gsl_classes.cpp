#include "gsl_classes.h"
/* #include <gsl/gsl_sf_result.h> */

GSLrng::GSLrng(unsigned long int seed)
{
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, seed);
}

GSLrng::~GSLrng()
{
	gsl_rng_free(r);
}
