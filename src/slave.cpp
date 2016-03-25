#include "mpi.h"
#include <math.h>
#include <gsl/gsl_const.h>
#include "ebeam.h"
#include "support_func.h"
#include "ions.h"
#include "fields.h"

int slave(int &p, int &id, MPI::Intracomm &slave_comm_id, bool verbose)
{
	int buf, step_buf;
	bool loop_alive;
	SimParams simparams_temp;

	/* const double m_e = GSL_CONST_MKSA_MASS_ELECTRON; */

	// ==================================
	// Receive starting gun
	// ==================================
	MPI::Status status;
	id = MPI::COMM_WORLD.Get_rank();

	MPI::COMM_WORLD.Recv(&buf, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
	if (verbose) printf("Slave %d says: **DUM DUM DUM DUM**\n", id);

	// ==================================
	// Receive run info
	// ==================================

	// ==================================
	// Make simparams
	// ==================================
	simparams_temp.bcast_receive();

	// ==================================
	// Recalculate to distribute
	// simulation across nodes
	// ==================================
	simparams_temp.n_e    /= (p-1);
	simparams_temp.n_ions /= (p-1);

	const SimParams simparams = simparams_temp;

	Field *field;
	field = new Field(simparams);

	// ==================================
	// Collectively create output file
	// ==================================
	ionsim::overwrite_file_parallel(simparams.filename, slave_comm_id);

	// ==================================
	// Write attributes
	// ==================================
	simparams.write_attributes_parallel(slave_comm_id);

	// ==================================
	// Generate beam
	// ==================================
	double nb_0, sr;

	Emit emit;
	emit.set_emit_n(simparams.emit_n, simparams.E);
	Plasma plas(simparams.n_p_cgs, simparams.m_ion_amu);
	Match mat(plas, simparams.E, emit);
	Beam x_beam(mat.beta(), mat.alpha(), emit);
	Beam y_beam(mat.beta(), mat.alpha(), emit);

	// if ((id == 1) && (verbose))
	// {
	// 	printf("Match- beta_x: %.5e\n", mat.beta());
	// 	printf("Sigma: %.5e\n", x_beam.sigma());
	// 	printf("Radius: %.5e\n", simparams.radius);
	// }

	double cov[2][2];
	x_beam.cov(cov);
	/* if (verbose) printf("Cov:\n[[ %.6e, %.6e ],\n [ %.6e, %.6e ]]\n", cov[0][0], cov[0][1], cov[1][0], cov[1][1]); */

	sr = x_beam.sigma();
	nb_0 = simparams.q_tot / (pow(2*M_PI, 1.5) * simparams.sz * sr * sr);

	Ebeam ebeam(simparams, x_beam, y_beam, id + 1);
	/* int lies   = 5; */
	/* int nolies = 1; */
	/* Ebeam ebeam(simparams, nolies, lies*lies, id + 1); */
	// Fix for having less charge per particle with more processors
	ebeam.qpp /= (p-1);

	if (id == 1) {
		printf("Qpp: %0.3e\n", ebeam.qpp);
		printf("Qtot: %0.3e\n", ebeam.qpp*ebeam.n_pts);
	}

	// ==================================
	// Generate ions
	// ==================================
	Ions ions(simparams, plas, simparams.n_ions, simparams.radius, simparams.length);

	// ==================================
	// Generate field
	// ==================================
	loop_alive = true;
	do
	{
		MPI::COMM_WORLD.Bcast(&buf, 1, MPI::INT, 0);
		switch (buf)
		{
			case ionsim::LOOP_KILL:
				loop_alive = false;
				break;

			case ionsim::LOOP_DUMP_IONS:
				MPI::COMM_WORLD.Bcast(&step_buf, 1, MPI::INT, 0);
				ions.dump_parallel(simparams.filename, step_buf, slave_comm_id);
				break;

			case ionsim::LOOP_DUMP_E:
				MPI::COMM_WORLD.Bcast(&step_buf, 1, MPI::INT, 0);
				ebeam.dump_parallel(simparams.filename, step_buf, slave_comm_id);
				break;

			case ionsim::LOOP_PUSH_IONS:
				switch (simparams.pushmethod)
				{
					case ionsim::PUSH_RUNGE_KUTTA:
						ions.push(simparams.dt, nb_0, sr);
						break;
					case ionsim::PUSH_SIMPLE:
						ions.push_simple(simparams.dt, nb_0, sr);
						break;
					case ionsim::PUSH_FIELD:
						delete field;
						field = new Field(simparams);

						(*field).recv_field(0);
						
						ions.push_field(simparams.dt, *field);

						break;
				}
				break;
			case ionsim::LOOP_GET_EFIELD:
				delete field;
				field = new Field(simparams);
				ebeam.field(*field);
				(*field).send_field(0);
				break;
		}

	} while ( loop_alive == true );


	return 0;
}
