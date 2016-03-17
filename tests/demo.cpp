#include "fields.h"
#include "consts.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
/* #include "support_func.h" */

double ** alloc_2d_array(long rowCount, long colCount);
int dealloc_2d_array(double &arr, long rowCount, long colCount);

int main(int argv, char **argc)
{
	printf("%f", sqrt(-1));
	long rowCount = 3;
	long colCount = 3;
	double ** arr;
	arr = alloc_2d_array(rowCount, colCount);
	arr[0][0] = 1;
	
	std::cout << "hi" << &arr << " then " << arr[0][0] << std::endl;
	dealloc_2d_array(**arr, rowCount, colCount);

	double ** out;
	out = new double*[rowCount];
	for (long i=0; i < rowCount; i++)
	{
		out[i] = new double[colCount];
	}
	out[0][0] = 1;

	/* out[1][2] = 3.14159; */
	/* out[1][2] *= 2; */
	/* std::cout << out[1][2] << std::endl; */

	/* for (long i=0; i < rowCount; i++) */
	/* { */
	/* 	delete [] out[i]; */
	/* } */
	/* delete [] out; */
	return 0;
}

double ** alloc_2d_array(long rowCount, long colCount)
{
	double lies = 1;
	double ** out = new double*[rowCount];
	for (long i=0; i < rowCount; i++)
	{
		out[i] = new double[colCount];
		printf("%ld\n", i);
	}
	out[0][0]=lies;

	std::cout << out[0][0] << std::endl;
	return out;
}

int dealloc_2d_array(double &arr, long rowCount, long colCount)
{
	/* long sizeY = sizeof(arr); */
	for (long i=0; i < rowCount; i++)
	{
		std::cout << i << std::endl;
		/* std::cout << arr[i][0] << std::endl; */
		/* delete [] arr[i]; */
	}
	/* delete [] arr; */
	return 0;
}

