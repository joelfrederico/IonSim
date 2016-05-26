#include "gtest/gtest.h" // we will add the path to C preprocessor later
#include "gtest.h"
#include "fieldtest.h"
#include "hdf5test.h"
#include "ionstest.h"
#include "simparamstest.h"

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
