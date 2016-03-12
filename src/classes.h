#ifndef __CLASSES_H_INCLUDED__
#define __CLASSES_H_INCLUDED__

#include "baseclass.h"

#include "mpi.h"
#include "beam.h"

/* class Ions; */
class Ebeam;

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
};

class Ions : public Parts
{
	private:
		double _radius;
		double _part_charge;

		Plasma * _plasma;

	public:
		Ions();
		Ions(Plasma * plasma, int n_pts, double radius, double length);

		int dump(std::string const &filename, int step, MPI::Intracomm &comm);
		int push(double dt, double nb_0, double sig_r);
		int push_simple(double dt, double nb_0, double sig_r);
};

#endif
