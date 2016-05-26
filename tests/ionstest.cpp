#include "ionstest.h"
#include "support.h"
#include "ions.h"

IonsTest::IonsTest()
{
	simparams = simparams_gen();
}

TEST_F(IonsTest, FieldForce)
{
	std::complex<double> Fr = F_r(2e-6, 2e-6, simparams);
	EXPECT_EQ(Fr.real(), 0);
	EXPECT_EQ(Fr.imag(), 0);
}

