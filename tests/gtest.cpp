#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "fields.h"

TEST_F(FieldTest, IsZeroInitially)
{
	long x_size = 10;
	long y_size = 10;
	complex_double test;
	bool pass = true;
	double re, im;

	Field field(x_size, y_size);
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			test = field(i, j);
			re = test.real();
			im = test.imag();
			if ( (re != 0) || (im != 0 ) )
			{
				pass = false;
				ADD_FAILURE() << "Complex number not zero at (" << i << ", " << j << "): (" << re << ", " << im << ")";
			}
		}
	}
	/* ASSERT_TRUE(pass) << "The field is not initialized to zero."; */
	/* ASSERT_EQ(field.get(0, 0).real(), 0) << "Here"; */
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
