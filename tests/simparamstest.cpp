#include "simparamstest.h"
#include "simparams.h"
#include <string.h>
#include "support.h"

TEST_F(SimparamsTest, dt)
{
	simparams = simparams_gen();
	EXPECT_NEAR(simparams.dt(), 6.34737e-14, 1e-17);
}

TEST_F(SimparamsTest, ion_mass)
{
	simparams = simparams_gen();
	EXPECT_NEAR(simparams.ion_mass(), 1.6737234599290802e-27, 1e-30);
}

TEST_F(SimparamsTest, LoadXML)
{
	std::string file = PKGDATADIR;
	simparams = SimParams(file + "/config.xml", false);
}
