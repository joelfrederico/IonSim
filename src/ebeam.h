#ifndef __EBEAM_H_INCLUDED__
#define __EBEAM_H_INCLUDED__

#include "parts.h"
#include "beam.h"
#include <mpi.h>
#include <complex>
#include "consts.h"
#include "field_data.h"

class Ebeam : public Parts
{
	private:
		// ==================================
		// Private member data
		// ==================================
		SimParams _simparams;

		int _gen_bivariate_gaussian(unsigned long int s, Cov x_cov, Cov y_cov, Cov z_cov);
		
		const Cov _x_cov;
		const Cov _y_cov;
		const Cov _z_cov;

		double n_resolve;
		double srsq;

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		Ebeam(const SimParams &simparams, Beam x_beam, Beam y_beam, unsigned long int s);
		Ebeam(const SimParams &simparams, double sx, double sxp, double sy, double syp, unsigned long int s);
		Ebeam(
				const SimParams &simparams,
				const double n_pts,
				const double type,
				double_vec x_in,
 				double_vec xp_in,
 				double_vec y_in,
 				double_vec yp_in,
 				double_vec z_in,
 				double_vec zp_in
				);

		// ==================================
		// Member data
		// ==================================
		double qpp;
		double z_end;

		// ==================================
		// Member methods
		// ==================================
		double x_std();
		double y_std();
		double x_mean();
		double y_mean();

		Ebeam between(double z0, double z1);

		int field_BE(Field_Data &field);
		int field_Coulomb(Field_Data &field);
		int field_Coulomb_sliced(Field_Data &field);
};

#endif
