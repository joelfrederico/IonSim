#ifndef __WRITER_SERIAL_H_INCLUDED__
#define __WRITER_SERIAL_H_INCLUDED__

#include "writer_base.h"
#include "parts.h"
#include "field_data.h"

class WriterSerial : public WriterBase
{
	private:
		int open_file(std::string const &filename);

	public:
		WriterSerial(const std::string &filename);

		int overwrite_file_serial(std::string const &filename);

		int writedata(long step, std::string const &group, std::string const &dataset, const Parts &parts);
		int writedata(long step, const Field_Data &field);
};

#endif
