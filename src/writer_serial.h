#ifndef __WRITER_SERIAL_H_INCLUDED__
#define __WRITER_SERIAL_H_INCLUDED__

#include "writer_base.h"
#include "parts.h"
#include "field_data.h"

class WriterSerial : public WriterBase
{
	private:
		int open_file();
		int overwrite_file_serial();

		int _init(const std::string &filename, bool overwrite);

	public:
		WriterSerial(const std::string &filename);
		WriterSerial(const std::string &filename, bool overwrite);


		int writedata(long step, const std::string &group, const std::string &dataset_str, const Parts &parts);
		int writedata(long step, Field_Data &field);
};

#endif
