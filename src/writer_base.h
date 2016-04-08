#ifndef __WRITER_BASE_H_INCLUDED__
#define __WRITER_BASE_H_INCLUDED__

#include <string>
#include <hdf5.h>
#include <typeinfo>
#include "field_data.h"
#include "simparams.h"
#include "hdf5_classes.h"

typedef int writer_t;
const writer_t WRITER_NULL     = 0;
const writer_t WRITER_SERIAL   = 1;
const writer_t WRITER_PARALLEL = 2;

class WriterBase
{
	protected:
		std::string _filename;
		hid_t file_id;
		writer_t writer_type;

		hid_t dataset_create(hid_t &group_id, hid_t &dataspace_id, int rank, hsize_t *count, std::string const &dataset);

	public:
		WriterBase(const std::string &filename);
		~WriterBase();

		int write_attributes(const SimParams &simparams) const;
};

#endif
