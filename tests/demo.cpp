#include "fields.h"
#include "consts.h"

int main(int argv, char **argc)
{
	long x_size = 3;
	long y_size = 5;
	long z_size = 7;
	Field field(x_size, y_size, z_size);
	for (int k=0; k < z_size; k++)
	{
		for (int j=0; j < y_size; j++)
		{
			for (int i=0; i < x_size; i++)
			{
				printf("%d, %d: (%0.3e, %0.3e)\n", i, j, field.x(i, j, k), field.x(i, j, k));
			}
		}
	}

	return 0;
}
