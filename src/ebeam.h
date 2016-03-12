#ifndef __EBEAM_H_INCLUDED__
#define __EBEAM_H_INCLUDED__

#include "baseclass.h"
#include "beam.h"
#include "mpi.h"

class Ebeam : public Parts
{
	private:
		double _q_tot;
		double _E;
		
		Beam _x_beam;
		Beam _y_beam;

		double x_cov[2][2];
		double y_cov[2][2];
		double z_cov[2][2];
	public:
		/// Create Ebeam class.
		Ebeam(int nparts, double mass, double q_tot, double E, Beam x_beam, Beam y_beam, double z_cov[2][2]);

		int dump(std::string const &filename, int step, MPI::Intracomm &comm);
		double x_std();
		double y_std();
};

#endif
