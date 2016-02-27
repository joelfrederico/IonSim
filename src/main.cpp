#include "config.h"
#include <stdio.h>

/* using namespace std; */

int main(int argc, char **argv)
{
	printf("Hello world\n");
	printf(PACKAGE_NAME);
	printf(" ");
	printf(PACKAGE_VERSION);
	return 0;
}
