#ifndef __SIMPARAMS_H_INCLUDED__
#define __SIMPARAMS_H_INCLUDED__

#include <iomanip>
#include "beam.h"
#include "emit.h"
#include <string>
#include <mpi.h>
#include "consts.h"
#include <typeinfo>
#include "loop_comm.h"
#include "pugixml/src/pugixml.hpp"
#include <stdexcept>

const int WIDTH = 20;

template <class T>
void getdata(pugi::xml_node node, std::string key, bool verbose, T &out)
{
	if (!node)
	{
		throw std::runtime_error("XML doesn't exist, suspect invalid config file. Key: " + key);
	} else {
		pugi::xml_text text = node.child(key.c_str()).text();
		if (typeid(out) == typeid(double)) {
			out = text.as_double();
		} else if (typeid(out) == typeid(long long)) {
			out = text.as_double();
		} else if (typeid(out) == typeid(int)) {
			out = text.as_int();
		} else {
			throw std::runtime_error("Type not matched");
		}
	}

	if (verbose) std::cout << std::left << std::setw(WIDTH) << key << ": " << out << std::endl;
}

class SimParams
{
	private:
		double *_dt;
		int _init();

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		SimParams(std::string filename);
		SimParams(std::string filename, bool verbose);
		SimParams();
		/* SimParams(const SimParams &obj); */
		~SimParams();

		// ==================================
		// Communications
		// ==================================
		int bcast_send() const;
		int bcast_receive();

		/* int write_attributes_parallel(Writer &writer) const; */

		template <class T>
		int bcast_send_wrap(T send) const
		{
			T buf = send;
			if (typeid(send) == typeid(double)) {
				MPI_Bcast(&buf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			} else if (typeid(send) == typeid(int)) {
				MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
			} else if (typeid(send) == typeid(long long)) {
				MPI_Bcast(&buf, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
			} else return -1;
			return 0;
		}

		// ==================================
		// Config
		// ==================================
		// - Generic
		std::string filename;

		// - Ion Window
		double field_trans_wind;
		double radius;
		double z_end;
		int n_field_x;
		int n_field_y;
		int n_field_z;

		// - Electron sim
		int n_steps;

		// - Beam
		double E;
		double emit_n;
		double q_tot;
		double sdelta;
		double sz;
		double z_center;
		long long n_e;
		zdist_t zdist;

		// - Ions
		double m_ion_amu;
		double n_p_cgs;
		long long n_ions;
		pushmethod_t pushmethod;

		// ==================================
		// Data methods
		// ==================================
		int z_cov(double (&out)[2][2]);
		double ion_mass() const;
		double dt() const;
		double gamma_rel() const;
		double dz() const;
		long long n_e_node() const;
		long long n_ions_node() const;
		double qpp_e() const;
};

#endif
