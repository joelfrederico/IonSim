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

		int _gen_bivariate_gaussian(unsigned long int s, double x_cov[2][2], double y_cov[2][2], double z_cov[2][2]);

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		Ebeam(const SimParams &simparams, Beam x_beam, Beam y_beam, unsigned long int s);
		Ebeam(const SimParams &simparams, double sx, double sy, unsigned long int s);
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

		// ==================================
		// Member methods
		// ==================================
		int dump_parallel(std::string const &filename, int step, MPI::Intracomm &comm);
		int dump_serial(std::string const &filename, int step);

		double x_std();
		double y_std();
		double x_mean();
		double y_mean();

		Ebeam between(double z0, double z1);
		int field(Field_Data &field);
};

#endif
