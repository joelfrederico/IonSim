#ifndef __SUPPORT_FUNC_H_INCLUDED__
#define __SUPPORT_FUNC_H_INCLUDED__

#include <mpi.h>
#include <string>
#include "parts.h"
#include <hdf5.h>
#include <gsl/gsl_const_mksa.h>
#include "fields.h"

namespace ionsim
{
	// ==================================
	// Methods
	// ==================================
	double GeV2gamma(double GeV);
	double gamma2GeV(double gamma);
	double gaussian(double mean, double sigma);
	hid_t open_file_parallel(std::string const &filename, MPI::Intracomm &slave_comm_id);
	hid_t open_file(std::string const &filename);
	double ** alloc_2d_array(long rowCount, long colCount);
	int dealloc_2d_array(double ** (&arr), long rowCount);
	int overwrite_file_parallel(std::string const &filename, MPI::Intracomm &slave_comm_id);
	int overwrite_file_serial(std::string const &filename);
	int sendloop(const int &message);
	int loop_get_fields(Field &field);
	int loop_push_ions(Field &field);
	int sendloop(const int &message, int step);
	hid_t group_access(hid_t &file_id, std::string const &group);
	hid_t group_access(hid_t &file_id, const char *group);

	template <class T>
	int writeattribute(hid_t file_id, std::string const &attr_name, T attr_value)
	{
		hid_t type_id = H5T_NATIVE_DOUBLE;
		const std::type_info &type = typeid(attr_value);
		if (typeid(attr_value) == typeid(double)) {
			type_id = H5T_NATIVE_DOUBLE;
		} else if (typeid(attr_value) == typeid(int)) {
			type_id = H5T_NATIVE_INT;
		} else if (typeid(attr_value) == typeid(long)) {
			type_id = H5T_NATIVE_LONG;
		} else {
			printf("Couldn't write attribute\n");
		}

		hid_t attr_id, dataspace_id;
		herr_t status;
		const hsize_t mysize = 1;

		// ==================================
		// Open the file
		// ==================================
		/* file_id = ionsim::open_file_parallel(filename, slave_comm_id); */

		// ==================================
		// Create an appropriate dataspace
		// ==================================
		dataspace_id = H5Screate(H5S_SIMPLE);
		status = H5Sset_extent_simple(dataspace_id, 1, &mysize, &mysize);

		// ==================================
		// Create and write attribute
		// ==================================
		attr_id = H5Acreate(file_id, attr_name.c_str(), type_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(dataspace_id);
		status = H5Awrite(attr_id, type_id, &attr_value);
		H5Aclose(attr_id);

		// ==================================
		// Clean up
		// ==================================
		/* H5Fclose(file_id); */

		return 0;
	}

	int dump_parallel(std::string const &filename, long step, std::string const &dataset, MPI::Intracomm &comm, const Parts &parts);
	int dump_serial(std::string const &filename, long step, std::string const &group, std::string const &dataset, const Parts &ebeam);
	int dump_serial(std::string const &filename, long step, std::string const &group, const Field &field);
	
	// ==================================
	// Consts
	// ==================================
	const double ELECTRON_REST_ENERGY = GSL_CONST_MKSA_MASS_ELECTRON * pow(GSL_CONST_MKSA_SPEED_OF_LIGHT, 2);

	const int TAG_LOOP_INSTRUCT = 100;
	const int TAG_FIELD         = 200;

	const loopflag_t LOOP_DUMP_E     = 3;
	const loopflag_t LOOP_DUMP_IONS  = 2;
	const loopflag_t LOOP_GET_EFIELD = 6;
	const loopflag_t LOOP_KILL       = 1;
	const loopflag_t LOOP_PUSH_E     = 4;
	const loopflag_t LOOP_PUSH_IONS  = 5;

	const parttype_t PARTS_E   = 2;
	const parttype_t PARTS_ION = 1;

	const pushmethod_t PUSH_RUNGE_KUTTA = 1;
	const pushmethod_t PUSH_SIMPLE      = 2;
	const pushmethod_t PUSH_FIELD       = 3;
}

#endif
