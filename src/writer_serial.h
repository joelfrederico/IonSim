#ifndef __WRITER_SERIAL_H_INCLUDED__
#define __WRITER_SERIAL_H_INCLUDED__

#include "writer_base.h"

class WriterSerial : public WriterBase
{
	private:
		int open_file(std::string const &filename);

	public:
		WriterSerial(const std::string &filename);

};

#endif
