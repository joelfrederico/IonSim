#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "fields.h"

FieldTest::FieldTest() : x_size(10), y_size(10), field(10, 10)
{
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			field(i, j) = complex_double(i*j, i*j*2);
		}
	}
}

bool expect_eq_complex(complex_double a, complex_double b)
{
	if ( (a.real() != b.real()) || (a.imag() != b.imag()) )
	{
		ADD_FAILURE() << "      Expected: (" << a.real() << ", " << a.imag() << ")" << std::endl << "To be equal to: (" << b.real() << ", " << b.imag() << ")";
		return false;
	} else {
		return true;
	}
}

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
			pass = expect_eq_complex(field(i, j), complex_double(0, 0));
		}
	}
}

TEST_F(FieldTest, AssigmentWorks)
{
	for (int i=0; i < x_size; i++)
	{
		for (int j=0; j < y_size; j++)
		{
			/* expect_eq_complex(field(i, j), complex_double(i*j, i*j*2)); */
		}
	}
	EXPECT_EQ(1, 2);
}



int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
