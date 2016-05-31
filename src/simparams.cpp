#include "simparams.h"
#include <gsl/gsl_const_mksa.h>
#include <math.h>
#include "support_func.h"
#include "consts.h"
#include "pugixml/src/pugixml.hpp"
#include <stdexcept>
#include <algorithm>

// ==================================
// XML Loading
// ==================================
std::string getstr(pugi::xml_node node, std::string name)
{
	std::string text;
	text = node.child(name.c_str()).text().as_string();

	text.erase(std::remove_if(text.begin(), text.end(), ::isspace), text.end());

	std::cout << std::left << std::setw(WIDTH) << name << ": " << text << std::endl;
	return text;
}

SimParams::SimParams(std::string xmlfile)
{
	pugi::xml_text text;
	std::string string;
	pugi::xml_document doc;

	pugi::xml_parse_result result = doc.load_file(xmlfile.c_str());

	if (!result) throw std::runtime_error("Couldn't load file: \"" + xmlfile+ "\"");

	pugi::xml_node beam = doc.child("config").child("Beam");
	pugi::xml_node ions = doc.child("config").child("Ions");

	getdata(beam, "E"         , E         ) ;
	getdata(beam, "emit_n"    , emit_n    ) ;
	getdata(beam, "length"    , length    ) ;
	getdata(beam, "m_ion_amu" , m_ion_amu ) ;
	getdata(beam, "n_p_cgs"   , n_p_cgs   ) ;
	getdata(beam, "q_tot"     , q_tot     ) ;
	getdata(beam, "radius"    , radius    ) ;
	getdata(beam, "sz"        , sz        ) ;
	getdata(beam, "sdelta"    , sdelta    ) ;

	string = getstr(beam, "zdist");
	if (string == "Gauss") {
		zdist = Z_DIST_GAUSS;
	} else if (string == "Flat") {
		zdist = Z_DIST_FLAT;
	} else {
		throw std::runtime_error ("Not a valid option for zdist:" + string);
	}

	getdata(beam, "n_steps", n_steps);

	string = getstr(beam, "pushmethod");
	if (string == "Simple") {
		pushmethod = PUSH_SIMPLE;
	} else if (string == "Field") {
		pushmethod = PUSH_FIELD;
	} else if (string == "RungeKutta") {
		pushmethod = PUSH_RUNGE_KUTTA;
	} else {
		throw std::runtime_error ("Not a valid option for zdist:" + string);
	}

	getdata(beam, "n_e", n_e);
	getdata(beam, "n_ions", n_ions);
	getdata(beam, "n_field_x", n_field_x);
	getdata(beam, "n_field_y", n_field_y);
	getdata(beam, "n_field_z", n_field_z);
	getdata(beam, "field_trans_wind", field_trans_wind);
	getdata(beam, "z_end", z_end);

	filename = getstr(beam, "filename");
}

// ==================================
// Constructors, Destructor
// ==================================
SimParams::SimParams(
	double _E,
	double _emit_n,
	double _length,
	double _m_ion_amu,
	double _n_p_cgs,
	double _q_tot,
	double _radius,
	double _sz,
	double _sdelta,
	zdist_t _zdist,
	int _n_steps,
	pushmethod_t _pushmethod,
	long long _n_e,
	int _n_field_x,
	int _n_field_y,
	int _n_field_z,
	double _field_trans_wind,
	double _z_end,
	long long _n_ions,
	std::string _filename
	)
{
	E                = _E;
	emit_n           = _emit_n;
	length           = _length;
	m_ion_amu        = _m_ion_amu;
	n_p_cgs          = _n_p_cgs;
	q_tot            = _q_tot;
	radius           = _radius;
	sz               = _sz;
	sdelta           = _sdelta;
	zdist            = _zdist;
	n_steps          = _n_steps;
	pushmethod       = _pushmethod;
	n_e              = _n_e;
	n_field_x        = _n_field_x;
	n_field_y        = _n_field_y;
	n_field_z        = _n_field_z;
	field_trans_wind = _field_trans_wind;
	z_end            = _z_end;
	n_ions           = _n_ions;
	filename         = _filename;

	_init();
}

SimParams::SimParams()
{
	_init();
}

int SimParams::_init()
{
	gamma_rel = ionsim::GeV2gamma(E);
	return 0;
}

SimParams::~SimParams()
{
	/* delete _dt; */
	/* std::cout << "Destructed _dt" << std::endl; */
}

int SimParams::bcast_send() const
{
	char *cbuf;
	long cbuf_l;

	// Send numerical parameters
	bcast_send_wrap(n_e              );
	bcast_send_wrap(n_ions           );
	bcast_send_wrap(q_tot            );
	bcast_send_wrap(radius           );
	bcast_send_wrap(length           );
	bcast_send_wrap(E                );
	bcast_send_wrap(emit_n           );
	bcast_send_wrap(n_p_cgs          );
	bcast_send_wrap(m_ion_amu        );
	bcast_send_wrap(sz               );
	bcast_send_wrap(sdelta           );
	bcast_send_wrap(zdist            );
	bcast_send_wrap(n_steps          );
	bcast_send_wrap(pushmethod       );
	bcast_send_wrap(n_field_x        );
	bcast_send_wrap(n_field_y        );
	bcast_send_wrap(n_field_z        );
	bcast_send_wrap(field_trans_wind );
	bcast_send_wrap(z_end            );

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
	MPI_Bcast(&n_e              , 1 , MPI_LONG_LONG , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_ions           , 1 , MPI_LONG_LONG , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&q_tot            , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&radius           , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&length           , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&E                , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&emit_n           , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_p_cgs          , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&m_ion_amu        , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&sz               , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&sdelta           , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&zdist            , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_steps          , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&pushmethod       , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_field_x        , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_field_y        , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_field_z        , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&field_trans_wind , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&z_end            , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);

	// Receive string

	MPI_Bcast(&cbuf_l, 1, MPI_LONG, 0, MPI_COMM_WORLD);

	cbuf = new char[cbuf_l];
	MPI_Bcast(cbuf, cbuf_l, MPI_CHAR, 0, MPI_COMM_WORLD);
	filename = std::string(cbuf);
	delete [] cbuf;

	_init();

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

double SimParams::dt() const
{
	return ( z_end / ((n_field_z-1) * GSL_CONST_MKSA_SPEED_OF_LIGHT) );
	return 0;
	if (*_dt == -1)
	{
		*_dt = ( z_end / ((n_field_z-1) * GSL_CONST_MKSA_SPEED_OF_LIGHT) );
	}

	return *_dt;
}
