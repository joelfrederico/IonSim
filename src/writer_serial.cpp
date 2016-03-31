#include "writer_serial.h"

WriterSerial::WriterSerial(const std::string &filename) : WriterBase(filename)
{
	open_file(filename);
}

int WriterSerial::open_file(std::string const &filename)
{
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	return 0;
}
