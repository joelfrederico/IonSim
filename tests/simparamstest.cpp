#include "simparamstest.h"
#include "simparams.h"
#include <string.h>
#include "support.h"

SimparamsTest::SimparamsTest()
{
	xmlfile = PKGDATADIR;
	xmlfile += "/config.xml";
}

TEST_F(SimparamsTest, dt)
{
	simparams = simparams_gen();
	EXPECT_NEAR(simparams.dt(), 6.2229127307165854e-14, 1e-17);
}

TEST_F(SimparamsTest, ion_mass)
{
	simparams = simparams_gen();
	EXPECT_NEAR(simparams.ion_mass(), 1.6737234599290802e-27, 1e-30);
}

TEST_F(SimparamsTest, LoadXML)
{
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(xmlfile.c_str());

	ASSERT_TRUE(result) << "File could not be loaded.";
}

TEST_F(SimparamsTest, ParseXML)
{
	simparams = SimParams(xmlfile, false);
}
