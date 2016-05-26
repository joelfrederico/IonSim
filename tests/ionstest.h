#ifndef __IONSTEST_H_INCLUDED__
#define __IONSTEST_H_INCLUDED__

#include "gtest/gtest.h"
#include "simparams.h"
#include "support.h"

class IonsTest : public :: testing::Test
{
	public:
		SimParams simparams;

		IonsTest();
};

#endif
