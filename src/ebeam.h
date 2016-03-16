#ifndef __EBEAM_H_INCLUDED__
#define __EBEAM_H_INCLUDED__

#include "baseclass.h"
#include "beam.h"
#include "mpi.h"
#include <complex>
#include "consts.h"
#include "fields.h"

class Ebeam : public Parts
{
	private:
		SimParams * _simparams;
		
		Beam _x_beam;
		Beam _y_beam;

		double x_cov[2][2];
		double y_cov[2][2];
		double z_cov[2][2];

		double i_to_x(long i);
		double j_to_y(long j);
		double k_to_z(long k);
	public:
		/// Create Ebeam class.
		Ebeam(SimParams &simparams, Beam x_beam, Beam y_beam);
		Ebeam(
				SimParams &simparams,
				double_vec x_in,
 				double_vec xp_in,
 				double_vec y_in,
 				double_vec yp_in,
 				double_vec z_in,
 				double_vec zp_in
				);

		double q_tot;

		int dump(std::string const &filename, int step, MPI::Intracomm &comm);

		double x_std();
		double y_std();
		double x_mean();
		double y_mean();

		Ebeam between(double z0, double z1);
		int get_field(Field &field);
};

#endif
