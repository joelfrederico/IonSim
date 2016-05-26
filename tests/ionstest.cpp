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
	EXPECT_NEAR(Fr.real(), -3.76203e-08, 1e-12);
	EXPECT_NEAR(Fr.imag(), -3.76203e-08, 1e-12);
}
