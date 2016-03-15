#include "fields.h"
#include "consts.h"

int main(int argv, char **argc)
{
	long x = 10;
	long y = 10;
	Field field(x, y);
	for (int i=0; i < 10; i++)
	{
		for (int j=0; j < 10; j++)
		{
			complex_double temp = field(i, j);
			printf("%d, %d: (%0.3e, %0.3e)\n", i, j, temp.real(), temp.imag());
		}
	}

	return 0;
}