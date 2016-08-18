#include "slave_loops.h"
#include "simparams.h"
#include "scalar_data.h"
#include "ebeam.h"
#include "scalar_data_comm.h"
#include "writer_parallel.h"
#include "ions.h"
#include "field_data.h"
#include "field_comm.h"

template<typename T>
int SL_get_rho(const unsigned int step_buf, const unsigned int substep_buf, const SimParams &simparams, ScalarData<T> &rho, const Ebeam &ebeam)
{
	// ==================================
	// Initialize variables
	// ==================================
	ScalarData_Comm scalarcomm;
	long double z0, z1;

	// ==================================
	// Histogram to find rho
	// ==================================
	z0 = step_buf * simparams.dz();
	z1 = (step_buf+1) * simparams.dz();
	ebeam.get_rho_dz(z0, z1, rho, simparams);

	// ==================================
	// Send rho to master
	// ==================================
	scalarcomm.send_scalar(rho, 0);

	return 0;
}

int SL_dump_ions(const SimParams &simparams, const LoopComm loopcomm, const unsigned int &step_buf, const unsigned int &substep_buf, const Ions &ions)
{
	// ==================================
	// Initialize variables
	// ==================================
	WriterParallel writer_p(simparams.filename, loopcomm.slave_comm);
	std::string subgroup;
	std::stringstream streamme;

	subgroup = "ions_steps";
	streamme.str("");
	streamme << "ions_" << std::setfill('0') << std::setw(4) << substep_buf;

	writer_p.writedata_substep(step_buf, substep_buf, streamme.str(), subgroup, ions);

	return 0;
}

int SL_dump_electrons(const SimParams &simparams, const LoopComm loopcomm, const unsigned int &step_buf, const Ebeam &ebeam)
{
	WriterParallel writer_p = WriterParallel(simparams.filename, loopcomm.slave_comm);
	writer_p.writedata(step_buf, "electrons", ebeam);

	return 0;
}

int SL_push_ions(const SimParams &simparams, const unsigned int &substep_buf, const Field_Data *field, Ions &ions)
{
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

	return 0;
}

int SL_get_efield(const SimParams &simparams, const Ebeam &ebeam, Field_Data *&field)
{
	Field_Comm fieldcomm;

	// ==================================
	// Send E-field to Master
	// ==================================
	delete field;
	field = new Field_Data(simparams);

	/* ebeam.field_Coulomb(*field); */
	ebeam.field_Coulomb_sliced(*field);
	fieldcomm.send_field(*field, 0);
				
	return 0;
}

int SL_send_efield(const SimParams &simparams, Field_Data *&field)
{
	Field_Comm fieldcomm;

	delete field;
	field = new Field_Data(simparams);
	fieldcomm.recv_field_copy(*field, 0);

	return 0;
}

int SL_get_ifield(const Field_Data *ion_field)
{
	Field_Comm fieldcomm;
	fieldcomm.send_field(*ion_field, 0);

	return 0;
}
