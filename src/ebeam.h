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

		int _gen_bivariate_gaussian(Cov x_cov, Cov y_cov, Cov z_cov);
		
		const Cov _x_cov;
		const Cov _y_cov;
		const Cov _z_cov;

		double srsq;

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		Ebeam(const SimParams &simparams, Beam x_beam, Beam y_beam);
		Ebeam(const SimParams &simparams, double sx, double sxp, double sy, double syp);
		Ebeam(
				const SimParams &simparams,
				const double n_pts,
				const double type,
				ldouble_vec x_in,
 				ldouble_vec xp_in,
 				ldouble_vec y_in,
 				ldouble_vec yp_in,
 				ldouble_vec z_in,
 				ldouble_vec zp_in
				);

		// ==================================
		// Member data
		// ==================================
		double n_resolve() const;
		double n_0() const;
		double sr_macro() const;

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

		bool verbose();
};

#endif
