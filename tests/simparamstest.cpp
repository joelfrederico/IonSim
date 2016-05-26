#include "simparamstest.h"
#include "simparams.h"
#include <string.h>
#include "support.h"

SimparamsTest::SimparamsTest()
{
	simparams = simparams_gen();
}

TEST_F(SimparamsTest, dt)
{
	EXPECT_NEAR(simparams.dt(), 6.34737e-14, 1e-17);
}

TEST_F(SimparamsTest, ion_mass)
{
	EXPECT_NEAR(simparams.ion_mass(), 1.6737234599290802e-27, 1e-30);
}
