#include "mpi.h"
#include <math.h>
#include <gsl/gsl_const.h>
#include "ebeam.h"
#include "support_func.h"
#include "ions.h"
#include "field_data.h"
#include "writer_serial.h"
#include "writer_parallel.h"
#include "loop_comm.h"

int slave(bool verbose)
{
	int buf, step_buf;
	bool loop_alive;
	SimParams simparams_temp;
	Field_Comm fieldcomm;
	LoopComm loopcomm;

	// ==================================
	// Receive starting gun
	// ==================================
	MPI_Status status;

	MPI_Recv(&buf, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if (verbose) printf("Slave %d says: **DUM DUM DUM DUM**\n", loopcomm.id);

	MPI_Barrier(MPI_COMM_WORLD);

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
	simparams_temp.n_e    /= (loopcomm.p-1);
	simparams_temp.n_ions /= (loopcomm.p-1);

	const SimParams simparams = simparams_temp;

	Field_Data *field;
	field = new Field_Data(simparams);

	// ==================================
	// Collectively create output file
	// ==================================
	/* ionsim::overwrite_file_parallel(simparams.filename, slave_comm_id); */

	// ==================================
	// Write attributes
	// ==================================
	WriterParallel *writer_p;
	writer_p = new WriterParallel(simparams.filename, loopcomm.slave_comm, true);
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

	Ebeam ebeam(simparams, x_beam, y_beam, loopcomm.id + 1);
	// Fix for having less charge per particle with more processors
	ebeam.qpp /= (loopcomm.p-1);

	// ==================================
	// Generate ions
	// ==================================
	Ions ions(simparams, plas, simparams.n_ions, simparams.radius, simparams.length);

	// ==================================
	// Slave loop
	// ==================================
	loop_alive = true;
	do
	{
		MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		switch (buf)
		{
			case ionsim::LOOP_KILL:
				// ==================================
				// Terminate Loop
				// ==================================
				loop_alive = false;
				break;

			case ionsim::LOOP_DUMP_IONS:
				// ==================================
				// Write ions to file
				// ==================================
				ionsim::sendloop(step_buf);
				writer_p = new WriterParallel(simparams.filename, loopcomm.slave_comm);

				/* (*writer_p).writedata_substep(step_buf, */ 

				delete writer_p;

				break;

			case ionsim::LOOP_DUMP_E:
				// ==================================
				// Write electrons to file
				// ==================================
				MPI_Bcast(&step_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				/* ebeam.dump_parallel(simparams.filename, step_buf, slave_comm_id); */
				break;

			case ionsim::LOOP_PUSH_IONS:
				// ==================================
				// Push ions
				// ==================================
				MPI_Bcast(&step_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				switch (simparams.pushmethod)
				{
					case ionsim::PUSH_RUNGE_KUTTA:
						ions.push(simparams.dt, nb_0, sr);
						break;
					case ionsim::PUSH_SIMPLE:
						ions.push_simple(simparams.dt, nb_0, sr);
						break;
					case ionsim::PUSH_FIELD:
						ions.push_field(simparams.dt, *field, step_buf);

						break;
				}
				break;
			case ionsim::LOOP_GET_EFIELD:
				// ==================================
				// Send E-field to Master
				// ==================================
				std::cout << "Getting field" << std::endl;

				delete field;
				field = new Field_Data(simparams);

				ebeam.field_Coulomb(*field);
				fieldcomm.send_field(*field, 0);

				break;
			case ionsim::LOOP_SEND_EFIELD:
				// ==================================
				// Retrieve field from Master
				// ==================================
				delete field;
				field = new Field_Data(simparams);
				fieldcomm.recv_field_copy(*field, 0);
				break;
		}

	} while ( loop_alive == true );

	return 0;
}
