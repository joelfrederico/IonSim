#ifndef __CLASSES_H_INCLUDED__
#define __CLASSES_H_INCLUDED__

#include "baseclass.h"
#include "consts.h"

#include "mpi.h"
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>

class Plasma;
class Beam;
class Ions;
class Ebeam;



class Beam
{
	private:
		double _alpha;
		double _beta;
		Emit _emit;
	public:
		Beam();
		Beam(double beta, double alpha, Emit emit);

		double alpha();
		double beta();
		void cov(double output[2][2]);

};



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
		Ebeam(int nparts, double q_tot, double E, Beam x_beam, Beam y_beam, double z_cov[2][2]);
		~Ebeam();

		int dump(std::string const &filename, MPI::Intracomm &comm);
};


class Plasma
{
	private:
		double _n_p;
		double _ion_mass_amu;
	public:
		Plasma();
		Plasma(double n_p_cgs, double ion_mass_amu);

		double n_p();
		double m();
		double w_p();
		double k_ion(double E);
};

class Match
{
	private:
		Plasma _plasma;
		double _E;
		Emit _emit;
	public:
		Match(Plasma plasma, double E, Emit emit);
		double beta();
		double alpha();
};

class Ions : public Parts
{
	private:
		double _x_window;
		double _y_window;

		double_vec _delta;

		Plasma _plasma;

	public:
		Ions();
		Ions(Plasma plasma, int n_pts, double x_window, double y_window);
		double_vec x();
		double_vec xp();
		double_vec y();
		double_vec yp();
		double_vec z();
		double_vec delta();

		int dump(std::string const &filename);
};

#endif
