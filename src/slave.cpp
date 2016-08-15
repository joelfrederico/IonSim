#include "consts.h"
#include "ebeam.h"
#include "emit.h"
#include "fftw_classes.h"
#include "field_data.h"
#include "ions.h"
#include "loop_comm.h"
#include "mpi.h"
#include "scalar_data_comm.h"
#include "support_func.h"
#include "writer_parallel.h"
#include "writer_serial.h"
#include <fftw3-mpi.h>
#include <gflags/gflags.h>
#include <gsl/gsl_const.h>
#include <iomanip>
#include <math.h>
#include <sstream>

DECLARE_bool(verbose);

int slave()
{
	// ==================================
	// FFTW variables
	// ==================================
	ptrdiff_t N0, N1;
	ScalarData_Comm scalarcomm;

	// ==================================
	// Other variables
	// ==================================
	std::stringstream streamme;
	std::string subgroup;
	int buf, step_buf, substep_buf;
	long double z0, z1;
	bool loop_alive;
	SimParams simparams_temp;
	Field_Comm fieldcomm;
	LoopComm loopcomm;

	// ==================================
	// Receive starting gun
	// ==================================
	MPI_Recv(&buf, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if (FLAGS_verbose) printf("Slave %d says: **DUM DUM DUM DUM**\n", loopcomm.id);

	MPI_Barrier(MPI_COMM_WORLD);


	// ==================================
	// Initialize FFTW
	// ==================================
	fftwl_mpi_init();
	fftwl_mpi_broadcast_wisdom(MPI_COMM_WORLD);

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

	// ==================================
	// Recalculate to distribute
	// simulation across nodes
	// ==================================
	ScalarData<ldouble> psi(simparams);
	ScalarData<ldouble> rho(simparams);
	ScalarData<std::complex<ldouble>> psi_k(simparams);

	Field_Data *field;
	Field_Data *ion_field;
	field = new Field_Data(simparams);
	ion_field = new Field_Data(simparams.n_field_x, simparams.n_field_y, 1, simparams.field_trans_wind, simparams.field_trans_wind, 0);

	// ==================================
	// Write attributes
	// ==================================
	WriterParallel *writer_p;

	// ==================================
	// Generate beam
	// ==================================
	/* double nb_0, sr; */

	Emit emit;
	emit.set_emit_n(simparams.emit_n, simparams.E);
	Plasma plas(simparams.n_p_cgs, simparams.m_ion_amu);
	Match mat(plas, simparams.E, emit);

	Beam x_beam(mat.beta(), mat.alpha(), emit);
	Beam y_beam(mat.beta(), mat.alpha(), emit);

	/* sr = ionsim::sr(simparams.emit_n, simparams.E, simparams.n_p_cgs, simparams.m_ion_amu); */
	/* nb_0 = ionsim::nb_0(simparams.q_tot, simparams.sz, sr); */

	Ebeam ebeam(simparams, x_beam, y_beam);
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
			case LOOP_START_E_ITER:
				loopcomm.recv_master(&step_buf);
				break;

			case LOOP_START_I_ITER:
				loopcomm.recv_master(&substep_buf);
				break;

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
				writer_p = new WriterParallel(simparams.filename, loopcomm.slave_comm);

				(*writer_p).writedata(step_buf, "electrons", ebeam);

				delete writer_p;

				break;

			case LOOP_PUSH_IONS:
				// ==================================
				// Push ions
				// ==================================
				switch (simparams.pushmethod)
				{
					case PUSH_RUNGE_KUTTA:
						/* ions.push(nb_0, sr); */
						break;
					case PUSH_SIMPLE:
						/* ions.push_simple(nb_0, sr); */
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

			case LOOP_GET_IFIELD:
				// ==================================
				// Retrieve current ion field
				// ==================================
				/* ions.field_Coulomb_sliced(*ion_field, substep_buf); */

				fieldcomm.send_field(*ion_field, 0);

				break;

			case LOOP_RESET_IFIELD:
				delete ion_field;
				ion_field = new Field_Data(simparams);
				break;

			case LOOP_GET_RHO:
				// ==================================
				// Histogram to find rho
				// ==================================
				z0 = step_buf * simparams.dz();
				z1 = (step_buf+1) * simparams.dz();
				ebeam.get_rho_dz(z0, z1, rho, simparams);

				// ==================================
				// Compute psi
				// ==================================
				psi = -rho;

				JTF_PRINTVAL(psi.ind(44, 78, 0));
				// ==================================
				// Send psi to master
				// ==================================
				scalarcomm.send_scalar(psi, 0);

				break;

			case LOOP_GET_FIELDS:
				// ==================================
				// Prepare for FFT
				// ==================================
				/* N0 = simparams.x_pts; */
				/* N1 = simparams.y_pts; */
				JTF_PRINT(Starting sub);
				psi_k = psifftw_base(simparams, loopcomm);

				JTF_PRINT(Heyi);

				break;
		}

	} while ( loop_alive == true );

	// ==================================
	// Clean up
	// ==================================
	delete field;
	fftwl_mpi_cleanup();


	return 0;
}
