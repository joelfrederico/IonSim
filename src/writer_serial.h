#ifndef __WRITER_SERIAL_H_INCLUDED__
#define __WRITER_SERIAL_H_INCLUDED__

#include "writer_base.h"
#include "parts.h"
#include "field_data.h"
#include <ionsim.h>

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

		template<typename T>
		int writedata(ScalarData<T> &buf, const std::string name) const
		{
			std::vector<unsigned long> size = {buf.x_pts, buf.y_pts};
			return writedata(buf.data, size, name);
		}

		template<typename T, typename T2>
		int writedata(const std::vector<std::complex<T>> &buf, const std::vector<T2> size, const std::string name) const
		{
			// ==================================
			// Initialize all variables
			// ==================================
			typename std::vector<std::complex<T>>::size_type bsize = buf.size();
			std::vector<typename std::vector<std::complex<T>>::size_type> _size = {size[0], size[1], 2};
			std::vector<double> dbuf;
			dbuf.resize(2*bsize);
			for (int i=0; i<bsize; i++)
			{
				dbuf[i] = buf[i].real();
				dbuf[i+bsize] = buf[i].imag();
			}

			return writedata(dbuf, _size, name);

		}

		template<typename T, typename T2>
		int writedata(const std::vector<T> &buf, const std::vector<T2> size, const std::string name) const
		{
			herr_t status;
			hid_t hdf5type;
			std::vector<double> dbuf;
			std::vector<double> fbuf;
			void *mpi_buf;
			// ==================================
			// Pair type with HDF5 type
			// ==================================
			if ((typeid(T) == typeid(long double)) || (typeid(T) == typeid(double)))
			{
				JTF_PRINT(Type is ldouble);
				// Set data type
				hdf5type = H5T_NATIVE_DOUBLE;

				// Copy data to properly-typed buffer
				dbuf.resize(buf.size());
				std::copy(buf.begin(), buf.end(), dbuf.begin());
				mpi_buf = dbuf.data();
			} else if (typeid(T) == typeid(float)) {
				hdf5type = H5T_NATIVE_FLOAT;
				fbuf.resize(buf.size());
				std::copy(buf.begin(), buf.end(), fbuf.begin());
				mpi_buf = fbuf.data();
			}

			// ==================================
			// Access or create a new group
			// ==================================
			/* GroupStepAccess step_group(file_id, step); */
			/* GroupAccess group(step_group.group_id, name); */
			/* GroupAccess group(file_id, name); */

			DatasetAccess vdataset(file_id, name, size, hdf5type);

			status = H5Dwrite(vdataset.dataset_id, hdf5type, vdataset.dataspace_id, vdataset.dataspace_id, H5P_DEFAULT, mpi_buf);

			return 0;
		}
};

#endif
