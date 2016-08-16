#ifndef __FFTW_CLASSES_H_INCLUDED__
#define __FFTW_CLASSES_H_INCLUDED__

#include "scalar_data.h"
#include <fftw3-mpi.h>

// ==================================
// Functions for sending data
// ==================================
int psifftw_base(LoopComm loopcomm);
int fftwl_recv_local_size(long long &local_n0, long long &local_0_start, const int id);

// ==================================
// Classes: receiving and doing FFT
// ==================================
class Fftw_base
{
	private:

	public:
		Fftw_base(const ptrdiff_t N0, const ptrdiff_t N1, LoopComm loopcomm);
};

class Fftw_2d_dft
{
	private:

	public:
		Fftw_2d_dft();

		int receive_data();
		int send_data();
		int execute();
};


class Fftw_2d_idft
{
	private:

	public:

};

#endif
