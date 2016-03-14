#ifndef __EBEAM_H_INCLUDED__
#define __EBEAM_H_INCLUDED__

#include "baseclass.h"
#include "beam.h"
#include "mpi.h"
#include <complex>
#include "consts.h"

class Ebeam : public Parts
{
	private:
		SimParams * _simparams;
		
		Beam _x_beam;
		Beam _y_beam;

		double x_cov[2][2];
		double y_cov[2][2];
		double z_cov[2][2];
	public:
		/// Create Ebeam class.
		Ebeam(SimParams &simparams, Beam x_beam, Beam y_beam);

		int dump(std::string const &filename, int step, MPI::Intracomm &comm);
		double x_std();
		double y_std();
		double x_mean();
		double y_mean();
};

#endif
