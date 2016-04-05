#include "mpi.h"
#include <math.h>
#include <gsl/gsl_const.h>
#include "ebeam.h"
#include "support_func.h"
#include "ions.h"
#include "field_data.h"
#include "writer_serial.h"
#include "writer_parallel.h"

int slave(int &p, int &id, MPI::Intracomm &slave_comm_id, bool verbose)
{
	int buf, step_buf;
	bool loop_alive;
	SimParams simparams_temp;
	Field_Comm fieldcomm;

	// ==================================
	// Receive starting gun
	// ==================================
	MPI::Status status;
	id = MPI::COMM_WORLD.Get_rank();

	MPI::COMM_WORLD.Recv(&buf, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
	if (verbose) printf("Slave %d says: **DUM DUM DUM DUM**\n", id);

	MPI::COMM_WORLD.Barrier();

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

	Field_Data *field;
	/* field = new Field_Data(simparams); */

	// ==================================
	// Collectively create output file
	// ==================================
	/* ionsim::overwrite_file_parallel(simparams.filename, slave_comm_id); */

	// ==================================
	// Write attributes
	// ==================================
	WriterParallel *writer_p;
	writer_p = new WriterParallel(simparams.filename, &slave_comm_id, true);
	(*writer_p).write_attributes(simparams);
	delete writer_p;

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

	double cov[2][2];
	x_beam.cov(cov);

	sr = x_beam.sigma();
	nb_0 = simparams.q_tot / (pow(2*M_PI, 1.5) * simparams.sz * sr * sr);

	Ebeam ebeam(simparams, x_beam, y_beam, id + 1);
	// Fix for having less charge per particle with more processors
	ebeam.qpp /= (p-1);

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
				/* ions.dump_parallel(simparams.filename, step_buf, slave_comm_id); */
				break;

			case ionsim::LOOP_DUMP_E:
				MPI::COMM_WORLD.Bcast(&step_buf, 1, MPI::INT, 0);
				/* ebeam.dump_parallel(simparams.filename, step_buf, slave_comm_id); */
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
						field = new Field_Data(simparams);

						/* (*field).recv_field(0); */
						
						ions.push_field(simparams.dt, *field);

						delete field;

						break;
				}
				break;
			case ionsim::LOOP_GET_EFIELD:
				std::cout << "Getting field" << std::endl;
				field = new Field_Data(simparams);
				ebeam.field_Coulomb(*field);
				fieldcomm.send_field(*field, 0);
				delete field;
				break;
		}

	} while ( loop_alive == true );

	return 0;
}
