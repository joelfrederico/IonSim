#include "mpi.h"
#include <math.h>
#include <gsl/gsl_const.h>
#include "ebeam.h"
#include "support_func.h"
#include "classes.h"

int slave(int &p, int &id, MPI::Intracomm &slave_comm_id)
{

	int buf;
	char *cbuf;
	long cbuf_l;
	int root = 0;

	long n_e;
	long n_ion;
	double q_tot;
	double radius;
	double length;
	double E;
	double emit_n;
	double n_p_cgs;
	double m_ion_amu;
	double sz;
	double sdelta;
	bool loop_alive = true;
	const double m_e = GSL_CONST_MKSA_MASS_ELECTRON;
	double t_tot;
	int n_steps;
	double dt;
	int runge_kutta;

	std::string filename;
	double nb_0, sr;

	// ==============================
	// Receive starting gun
	// ==============================
	MPI::Status status;

	MPI::COMM_WORLD.Recv(&buf, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
	printf("Slave %d says: **DUM DUM DUM DUM**\n", id);

	// ==============================
	// Receive run info
	// ==============================
	MPI::COMM_WORLD.Bcast(&n_e         , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&n_ion       , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&q_tot       , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&radius      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&length      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&E           , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&emit_n      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&n_p_cgs     , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&m_ion_amu   , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&sz          , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&sdelta      , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&t_tot       , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&n_steps     , 1 , MPI::INT    , 0);
	MPI::COMM_WORLD.Bcast(&dt          , 1 , MPI::DOUBLE , 0);
	MPI::COMM_WORLD.Bcast(&runge_kutta , 1 , MPI::INT    , 0);

	MPI::COMM_WORLD.Bcast(&cbuf_l, 1, MPI::LONG, 0);

	cbuf = new char[cbuf_l];
	MPI::COMM_WORLD.Bcast(cbuf, cbuf_l, MPI::CHAR, 0);
	filename = std::string(cbuf);
	delete [] cbuf;

	// ==============================
	// Create output file
	// ==============================
	ionsim::overwrite_file(filename, slave_comm_id);

	// ==============================
	// Write attributes
	// ==============================
	ionsim::writeattribute("n_e"       , n_e       , filename , slave_comm_id);
	ionsim::writeattribute("n_ion"     , n_ion     , filename , slave_comm_id);
	ionsim::writeattribute("q_tot"     , q_tot     , filename , slave_comm_id);
	ionsim::writeattribute("radius"    , radius    , filename , slave_comm_id);
	ionsim::writeattribute("length"    , length    , filename , slave_comm_id);
	ionsim::writeattribute("E"         , E         , filename , slave_comm_id);
	ionsim::writeattribute("emit_n"    , emit_n    , filename , slave_comm_id);
	ionsim::writeattribute("n_p_cgs"   , n_p_cgs   , filename , slave_comm_id);
	ionsim::writeattribute("m_ion_amu" , m_ion_amu , filename , slave_comm_id);
	ionsim::writeattribute("sz"        , sz        , filename , slave_comm_id);
	ionsim::writeattribute("sdelta"    , sdelta    , filename , slave_comm_id);
	ionsim::writeattribute("dt"        , dt        , filename , slave_comm_id);
	ionsim::writeattribute("n_steps"   , n_steps   , filename , slave_comm_id);

	// ==============================
	// Recalculate to distribute
	// simulation across nodes
	// ==============================
	n_e   /= (p-1);
	n_ion /= (p-1);

	// ==============================
	// Generate beam
	// ==============================
	double z_cov[2][2] = {{pow(sz, 2), 0}, {0, pow(sdelta, 2)}};

	Emit emit;
	emit.set_emit_n(emit_n, E);
	Plasma plas(n_p_cgs, m_ion_amu);
	Match mat(plas, E, emit);
	Beam x_beam(mat.beta(), mat.alpha(), emit);
	Beam y_beam(mat.beta(), mat.alpha(), emit);

	if (id == 1)
	{
		printf("Match- beta_x: %.5e\n", mat.beta());
		printf("Sigma: %.5e\n", x_beam.sigma());
		printf("Radius: %.5e\n", radius);
	}

	double cov[2][2];
	x_beam.cov(cov);
	printf("Cov:\n[[ %.6e, %.6e ],\n [ %.6e, %.6e ]]\n", cov[0][0], cov[0][1], cov[1][0], cov[1][1]);

	sr = x_beam.sigma();
	nb_0 = q_tot / (pow(2*M_PI, 1.5) * sz * sr * sr);

	Ebeam ebeam(n_e, m_e, q_tot, E, x_beam, y_beam, z_cov);

	// ==============================
	// Generate ions
	// ==============================
	Ions ions(&plas, n_ion, radius, length);

	do
	{
		MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
		switch (buf)
		{
			case ionsim::LOOP_KILL:
				loop_alive = false;
				break;

			case ionsim::LOOP_DUMP_IONS:
				MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
				ions.dump(filename, buf, slave_comm_id);
				break;

			case ionsim::LOOP_DUMP_E:
				MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
				ebeam.dump(filename, buf, slave_comm_id);
				break;

			case ionsim::LOOP_PUSH_IONS:
				if (runge_kutta == 1) {
					ions.push(dt, nb_0, sr);
				} else {
					ions.push_simple(dt, nb_0, sr);
				}
				break;
		}

	} while ( loop_alive == true );


	return 0;
}
