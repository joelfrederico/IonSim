#include "fields.h"
#include "consts.h"

int main(int argv, char **argc)
{
	const long X_PTS = 3;
	const long Y_PTS = 5;
	const long Z_PTS = 7;
	
	const double X_EDGE_MAG = 30e-6;
	const double Y_EDGE_MAG = 30e-6;
	const double Z_EDGE_MAG = 30e-6;

	Field field(X_PTS, Y_PTS, Z_PTS, X_EDGE_MAG, Y_EDGE_MAG, Z_EDGE_MAG);
	for (int k=0; k < Z_PTS; k++)
	{
		for (int j=0; j < Y_PTS; j++)
		{
			for (int i=0; i < X_PTS; i++)
			{
				printf("%d, %d: (%0.3e, %0.3e)\n", i, j, field.x(i, j, k), field.x(i, j, k));
			}
		}
	}

	return 0;
}
