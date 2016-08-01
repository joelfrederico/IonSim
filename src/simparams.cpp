#include "consts.h"
#include "pugixml/src/pugixml.hpp"
#include "simparams.h"
#include "support_func.h"
#include <algorithm>
#include <gsl/gsl_const_mksa.h>
#include <math.h>
#include <stdexcept>

// ==================================
// XML Loading
// ==================================
std::string getstr(pugi::xml_node node, std::string name, bool verbose)
{
	std::string text;
	text = node.child(name.c_str()).text().as_string();

	text.erase(std::remove_if(text.begin(), text.end(), ::isspace), text.end());

	if (verbose) std::cout << std::left << std::setw(WIDTH) << name << ": " << text << std::endl;

	return text;
}

int _register_file(SimParams &simparams, std::string xmlfile, bool verbose)
{
	pugi::xml_text text;
	std::string string;
	pugi::xml_document doc;

	pugi::xml_parse_result result = doc.load_file(xmlfile.c_str());

	if (!result) throw std::runtime_error("Couldn't load file: \"" + xmlfile+ "\"");

	if (verbose)
	{
		std::cout << "=========================================" << std::endl;
		std::cout << "Config file loaded" << std::endl;
		std::cout << "-----------------------------------------" << std::endl;
	}

	pugi::xml_node generic     = doc.child("config").child("Generic");
	pugi::xml_node ionsim      = doc.child("config").child("IonSim");
	pugi::xml_node electronsim = doc.child("config").child("ElectronSim");
	pugi::xml_node beam        = doc.child("config").child("Beam");
	pugi::xml_node ions        = doc.child("config").child("Ions");

	// ==================================
	// Generic
	// ==================================
	simparams.filename = getstr(generic, "filename", verbose);
	getdata(ionsim , "push_electrons" , verbose , simparams.push_electrons );
	getdata(ionsim , "push_ions"      , verbose , simparams.push_ions      );

	// ==================================
	// IonSim
	// ==================================
	getdata(ionsim , "field_trans_wind" , verbose , simparams.field_trans_wind );
	getdata(ionsim , "radius"           , verbose , simparams.radius           );
	getdata(ionsim , "n_field_x"        , verbose , simparams.n_field_x        );
	getdata(ionsim , "n_field_y"        , verbose , simparams.n_field_y        );
	getdata(ionsim , "n_field_z"        , verbose , simparams.n_field_z        );
	getdata(ionsim , "z_end"            , verbose , simparams.z_end            );
	getdata(ionsim , "ion_z_bool"       , verbose , simparams.ion_z_bool       );

	// ==================================
	// ElectronSim
	// ==================================
	getdata(electronsim, "n_steps", verbose, simparams.n_steps);
	
	// ==================================
	// Beam
	// ==================================
	getdata(beam , "n_e"      , verbose , simparams.n_e      );
	getdata(beam , "E"        , verbose , simparams.E        );
	getdata(beam , "emit_n"   , verbose , simparams.emit_n   );
	getdata(beam , "q_tot"    , verbose , simparams.q_tot    );
	getdata(beam , "sz"       , verbose , simparams.sz       );
	getdata(beam , "z_center" , verbose , simparams.z_center );
	getdata(beam , "sdelta"   , verbose , simparams.sdelta   );

	string = getstr(beam, "zdist", verbose);
	if (string == "Gauss") {
		simparams.zdist = Z_DIST_GAUSS;
	} else if (string == "Flat") {
		simparams.zdist = Z_DIST_FLAT;
	} else {
		throw std::runtime_error ("Not a valid option for zdist:" + string);
	}

	// ==================================
	// Ions
	// ==================================
	getdata(ions , "n_ions"    , verbose , simparams.n_ions    );
	getdata(ions , "m_ion_amu" , verbose , simparams.m_ion_amu );
	getdata(ions , "n_p_cgs"   , verbose , simparams.n_p_cgs   );

	string = getstr(ions, "pushmethod", verbose);
	if (string == "Simple") {
		simparams.pushmethod = PUSH_SIMPLE;
	} else if (string == "Field") {
		simparams.pushmethod = PUSH_FIELD;
	} else if (string == "RungeKutta") {
		simparams.pushmethod = PUSH_RUNGE_KUTTA;
	} else {
		throw std::runtime_error ("Not a valid option for zdist:" + string);
	}

	if (verbose)
	{
		std::cout << "=========================================" << std::endl;
	}

	return 0;
}

SimParams::SimParams(std::string xmlfile, bool verbose)
{
	_register_file(*this, xmlfile, verbose);
}

SimParams::SimParams(std::string xmlfile)
{
	_register_file(*this, xmlfile, true);
}

// ==================================
// Constructors, Destructor
// ==================================
SimParams::SimParams()
{
}

int SimParams::_init()
{
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
	bcast_send_wrap(E                );
	bcast_send_wrap(emit_n           );
	bcast_send_wrap(n_p_cgs          );
	bcast_send_wrap(m_ion_amu        );
	bcast_send_wrap(sz               );
	bcast_send_wrap(z_center         );
	bcast_send_wrap(sdelta           );
	bcast_send_wrap(zdist            );
	bcast_send_wrap(n_steps          );
	bcast_send_wrap(pushmethod       );
	bcast_send_wrap(n_field_x        );
	bcast_send_wrap(n_field_y        );
	bcast_send_wrap(n_field_z        );
	bcast_send_wrap(field_trans_wind );
	bcast_send_wrap(z_end            );
	bcast_send_wrap(push_electrons   );
	bcast_send_wrap(push_ions        );
	bcast_send_wrap(ion_z_bool       );

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
	MPI_Bcast(&E                , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&emit_n           , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_p_cgs          , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&m_ion_amu        , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&sz               , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&z_center         , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&sdelta           , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&zdist            , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_steps          , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&pushmethod       , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_field_x        , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_field_y        , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&n_field_z        , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&field_trans_wind , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&z_end            , 1 , MPI_DOUBLE    , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&push_electrons   , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&push_ions        , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);
	MPI_Bcast(&ion_z_bool       , 1 , MPI_INT       , 0 , MPI_COMM_WORLD);

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
	return dz() / GSL_CONST_MKSA_SPEED_OF_LIGHT;
}

double SimParams::dz() const
{
	return z_end / n_field_z;
}

double SimParams::gamma_rel() const
{
	return ionsim::GeV2gamma(E);
}

long long SimParams::n_e_node() const
{
	int p;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	return n_e / (p-1);
}

long long SimParams::n_ions_node() const
{
	int p;
	double out;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	out = n_ions / (p-1);
	return out;
}

double SimParams::qpp_e() const
{
	return q_tot / n_e;
}
