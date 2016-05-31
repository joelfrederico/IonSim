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
void getdata(pugi::xml_node node, std::string key, T &out)
{
	if (!node)
	{
		throw std::runtime_error("Xml doesn't exist, suspect invalid config file");
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

	std::cout << std::left << std::setw(WIDTH) << key << ": " << out << std::endl;
}

class SimParams
{
	private:
		int _init();
		double *_dt;

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		SimParams(std::string filename);
		SimParams(
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
			);

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
		// Data
		// ==================================
		double E;
		double emit_n;
		double field_trans_wind;
		double length;
		double m_ion_amu;
		double n_p_cgs;
		double q_tot;
		double radius;
		double sdelta;
		double sz;
		double z_end;
		int n_field_x;
		int n_field_y;
		int n_field_z;
		int n_steps;
		long long n_e;
		long long n_ions;

		pushmethod_t pushmethod;
		std::string filename;
		zdist_t zdist;

		// ==================================
		// Calculated Data
		// ==================================
		double gamma_rel;

		// ==================================
		// Data methods
		// ==================================
		int z_cov(double (&out)[2][2]);
		double ion_mass() const;
		double dt() const;
};

#endif
