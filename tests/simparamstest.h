#ifndef __SIMPARAMSTEST_H_INCLUDED__
#define __SIMPARAMSTEST_H_INCLUDED__

#include "gtest/gtest.h"
#include "simparams.h"

class SimparamsTest : public :: testing::Test
{
	private:

	public:
		SimparamsTest();

		SimParams simparams;
		std::string xmlfile;
};

#endif
