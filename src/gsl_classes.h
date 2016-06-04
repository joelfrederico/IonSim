#ifndef __GSL_CLASSES_H_INCLUDED__
#define __GSL_CLASSES_H_INCLUDED__

#include <gsl/gsl_rng.h>
/* #include <gsl/gsl_sf_result.h> */

class GSLrng
{
	private:
		
	public:
		GSLrng(unsigned long int seed);
		~GSLrng();

		gsl_rng *r;
};

#endif
