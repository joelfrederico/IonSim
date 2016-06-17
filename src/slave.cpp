#include "mpi.h"
#include "emit.h"
#include <math.h>
#include <gsl/gsl_const.h>
#include "ebeam.h"
#include "support_func.h"
#include "ions.h"
#include "field_data.h"
#include "writer_serial.h"
#include "writer_parallel.h"
#include "loop_comm.h"
#include "consts.h"
#include <sstream>
#include <iomanip>
#include <gflags/gflags.h>

DECLARE_bool(verbose);

int slave()
{
	std::stringstream streamme;
	std::string subgroup;
	int buf, step_buf, substep_buf;
	bool loop_alive;
	SimParams simparams_temp;
	Field_Comm fieldcomm;
	LoopComm loopcomm;

	// ==================================
	// Receive starting gun
	// ==================================
	MPI_Status status;

	MPI_Recv(&buf, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if (FLAGS_verbose) printf("Slave %d says: **DUM DUM DUM DUM**\n", loopcomm.id);

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
	/* simparams_temp.n_e    /= (loopcomm.p-1); */
	/* simparams_temp.n_ions /= (loopcomm.p-1); */

	const SimParams simparams = simparams_temp;

	Field_Data *field;
	field = new Field_Data(simparams);

	// ==================================
	// Write attributes
	// ==================================
	WriterParallel *writer_p;

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

	sr = ionsim::sr(simparams.emit_n, simparams.E, simparams.n_p_cgs, simparams.m_ion_amu);
	nb_0 = ionsim::nb_0(simparams.q_tot, simparams.sz, sr);

	Ebeam ebeam(simparams, x_beam, y_beam, loopcomm.id + 1);
	// Fix for having less charge per particle with more processors
	/* ebeam.qpp /= (loopcomm.p-1); */

	// ==================================
	// Generate ions
	// ==================================
	Ions ions(&simparams, plas);

	// ==================================
	// Slave loop
	// ==================================
	loop_alive = true;
	do
	{
		loopcomm.instruct(&buf);
		switch (buf)
		{
			case LOOP_KILL:
				// ==================================
				// Terminate Loop
				// ==================================
				loop_alive = false;
				break;

			case LOOP_DUMP_IONS:
				// ==================================
				// Write ions to file
				// ==================================
				loopcomm.recv_master(&step_buf);
				loopcomm.recv_master(&substep_buf);
				writer_p = new WriterParallel(simparams.filename, loopcomm.slave_comm);

				subgroup = "ions_steps";
				streamme.str("");
				streamme << "ions_" << std::setfill('0') << std::setw(4) << substep_buf;

				(*writer_p).writedata_substep(step_buf, substep_buf, streamme.str(), subgroup, ions);

				delete writer_p;

				break;

			case LOOP_DUMP_E:
				// ==================================
				// Write electrons to file
				// ==================================
				loopcomm.recv_master(&step_buf);

				writer_p = new WriterParallel(simparams.filename, loopcomm.slave_comm);

				(*writer_p).writedata(step_buf, "electrons", ebeam);

				delete writer_p;

				break;

			case LOOP_PUSH_IONS:
				// ==================================
				// Push ions
				// ==================================
				loopcomm.recv_master(&substep_buf);
				switch (simparams.pushmethod)
				{
					case PUSH_RUNGE_KUTTA:
						ions.push(nb_0, sr);
						break;
					case PUSH_SIMPLE:
						ions.push_simple(nb_0, sr);
						break;
					case PUSH_FIELD:
						ions.push_field(*field, substep_buf);

						break;
				}
				break;
			case LOOP_GET_EFIELD:
				// ==================================
				// Send E-field to Master
				// ==================================
				delete field;
				field = new Field_Data(simparams);

				/* ebeam.field_Coulomb(*field); */
				ebeam.field_Coulomb_sliced(*field);
				fieldcomm.send_field(*field, 0);

				break;
			case LOOP_SEND_EFIELD:
				// ==================================
				// Retrieve field from Master
				// ==================================
				delete field;
				field = new Field_Data(simparams);
				fieldcomm.recv_field_copy(*field, 0);
				break;
		}

	} while ( loop_alive == true );

	delete field;

	return 0;
}
