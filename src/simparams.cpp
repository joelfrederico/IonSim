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

	gamma_rel   = ionsim::GeV2gamma(_E);

	_init();
}

SimParams::SimParams()
{
	_init();
}

int SimParams::_init()
{
	return 0;
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
	MPI_Bcast(&cbuf_l, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	cbuf = new char[cbuf_l];
	strcpy(cbuf, filename.c_str());
	MPI_Bcast(cbuf, cbuf_l, MPI_CHAR, 0, MPI_COMM_WORLD);
	delete [] cbuf;

	return 0;
}

int SimParams::bcast_receive()
{
	long cbuf_l;
	char *cbuf;

	// Receive numerical parameters
	MPI_Bcast(&n_e        , 1 , MPI_LONG   , 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_ions     , 1 , MPI_LONG   , 0, MPI_COMM_WORLD);
	MPI_Bcast(&q_tot      , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&radius     , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&length     , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&E          , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&emit_n     , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_p_cgs    , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&m_ion_amu  , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&sz         , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&sdelta     , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&t_tot      , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_steps    , 1 , MPI_INT    , 0, MPI_COMM_WORLD);
	MPI_Bcast(&dt         , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD);
	MPI_Bcast(&pushmethod , 1 , MPI_INT    , 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_field_x  , 1 , MPI_LONG   , 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_field_y  , 1 , MPI_LONG   , 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_field_z  , 1 , MPI_LONG   , 0, MPI_COMM_WORLD);

	// Receive string

	MPI_Bcast(&cbuf_l, 1, MPI_LONG, 0, MPI_COMM_WORLD);

	cbuf = new char[cbuf_l];
	MPI_Bcast(cbuf, cbuf_l, MPI_CHAR, 0, MPI_COMM_WORLD);
	filename = std::string(cbuf);
	delete [] cbuf;

	gamma_rel   = ionsim::GeV2gamma(E);

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
