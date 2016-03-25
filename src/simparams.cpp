#include "simparams.h"
#include <gsl/gsl_const_mksa.h>
#include <math.h>
#include "support_func.h"
#include "consts.h"

// ==================================
// Constructors, Destructor
// ==================================
SimParams::SimParams(
	double _E,
	double _dt,
	double _emit_n,
	double _length,
	double _m_ion_amu,
	double _n_p_cgs,
	double _q_tot,
	double _radius,
	double _sdelta,
	double _sz,
	double _t_tot,
	int _n_steps,
	pushmethod_t _pushmethod,
	long _n_e,
	long _n_field_x,
	long _n_field_y,
	long _n_field_z,
	long _n_ions,
	std::string _filename
	)
{
	E           = _E;
	dt          = _dt;
	emit_n      = _emit_n;
	length      = _length;
	m_ion_amu   = _m_ion_amu;
	n_p_cgs     = _n_p_cgs;
	q_tot       = _q_tot;
	radius      = _radius;
	sdelta      = _sdelta;
	sz          = _sz;
	t_tot       = _t_tot;
	n_steps     = _n_steps;
	pushmethod = _pushmethod;
	n_e         = _n_e;
	n_field_x   = _n_field_x;
	n_field_y   = _n_field_y;
	n_field_z   = _n_field_z;
	n_ions      = _n_ions;
	filename    = _filename;
}

SimParams::SimParams()
{
}

int SimParams::bcast_send() const
{
	char *cbuf;
	long cbuf_l;

	// Send numerical parameters
	bcast_send_wrap(n_e        );
	bcast_send_wrap(n_ions     );
	bcast_send_wrap(q_tot      );
	bcast_send_wrap(radius     );
	bcast_send_wrap(length     );
	bcast_send_wrap(E          );
	bcast_send_wrap(emit_n     );
	bcast_send_wrap(n_p_cgs    );
	bcast_send_wrap(m_ion_amu  );
	bcast_send_wrap(sz         );
	bcast_send_wrap(sdelta     );
	bcast_send_wrap(t_tot      );
	bcast_send_wrap(n_steps    );
	bcast_send_wrap(dt         );
	bcast_send_wrap(pushmethod);
	bcast_send_wrap(n_field_x  );
	bcast_send_wrap(n_field_y  );
	bcast_send_wrap(n_field_z  );

	// Send string

	cbuf_l = filename.length()+1;
	MPI::COMM_WORLD.Bcast(&cbuf_l, 1, MPI::LONG, 0);
	cbuf = new char[cbuf_l];
	strcpy(cbuf, filename.c_str());
	MPI::COMM_WORLD.Bcast(cbuf, cbuf_l, MPI::CHAR, 0);
	delete [] cbuf;

	return 0;
}

int SimParams::bcast_receive()
{
	long cbuf_l;
	char *cbuf;

	// Receive numerical parameters
	MPI::COMM_WORLD.Bcast(&n_e         , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&n_ions      , 1 , MPI::LONG   , 0);
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
	MPI::COMM_WORLD.Bcast(&pushmethod , 1 , MPI::INT    , 0);
	MPI::COMM_WORLD.Bcast(&n_field_x   , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&n_field_y   , 1 , MPI::LONG   , 0);
	MPI::COMM_WORLD.Bcast(&n_field_z   , 1 , MPI::LONG   , 0);

	// Receive string

	MPI::COMM_WORLD.Bcast(&cbuf_l, 1, MPI::LONG, 0);

	cbuf = new char[cbuf_l];
	MPI::COMM_WORLD.Bcast(cbuf, cbuf_l, MPI::CHAR, 0);
	filename = std::string(cbuf);
	delete [] cbuf;
	return 0;
}

int SimParams::write_attributes_parallel(MPI::Intracomm &slave_comm_id) const
{
	hid_t file_id;
	file_id = ionsim::open_file_parallel(filename, slave_comm_id);

	ionsim::writeattribute(file_id, "n_e"       , n_e      );
	ionsim::writeattribute(file_id, "n_ions"    , n_ions   );
	ionsim::writeattribute(file_id, "q_tot"     , q_tot    );
	ionsim::writeattribute(file_id, "radius"    , radius   );
	ionsim::writeattribute(file_id, "length"    , length   );
	ionsim::writeattribute(file_id, "E"         , E        );
	ionsim::writeattribute(file_id, "emit_n"    , emit_n   );
	ionsim::writeattribute(file_id, "n_p_cgs"   , n_p_cgs  );
	ionsim::writeattribute(file_id, "m_ion_amu" , m_ion_amu);
	ionsim::writeattribute(file_id, "sz"        , sz       );
	ionsim::writeattribute(file_id, "sdelta"    , sdelta   );
	ionsim::writeattribute(file_id, "dt"        , dt       );
	ionsim::writeattribute(file_id, "n_steps"   , n_steps  );
	ionsim::writeattribute(file_id, "n_field_x" , n_field_x);
	ionsim::writeattribute(file_id, "n_field_y" , n_field_y);
	ionsim::writeattribute(file_id, "n_field_z" , n_field_z);
	H5Fclose(file_id);

	return 0;
}
// ==================================
// Public methods
// ==================================
int SimParams::z_cov(double (&out)[2][2])
{
	out[0][0] = pow(sz, 2);
	out[0][1] = out[1][0] = 0;
	out[1][1] = pow(sdelta, 2);
	return 0;
}

double SimParams::ion_mass() const
{
	return m_ion_amu * GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
}
