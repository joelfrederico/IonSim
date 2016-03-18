#include "ebeam.h"
#include "fields.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "support_func.h"

int main(int argv, char **argc)
{
	std::string filename = "test.h5";
	ionsim::overwrite_file_serial(filename);
	return 0;
}
