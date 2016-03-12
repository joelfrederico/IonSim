#include "beam.h"
#include <math.h>

// ==============================
// Beam
// ==============================
Beam::Beam() {}

Beam::Beam(double beta, double alpha, Emit emit)
{
	_beta  = beta;
	_alpha = alpha;
	_emit = emit;
}

double Beam::alpha()
{
	return _alpha;
}

double Beam::beta()
{
	return _beta;
}

double Beam::sigma()
{
	return sqrt(_beta*_emit.emit());
}

void Beam::cov(double output[2][2])
{
	output[0][0] = beta();
	output[0][1] = output[1][0] = -alpha();
	output[1][1] = (1.0+pow(alpha(), 2))/beta();
	double buf = _emit.emit();
	for (int i=0; i < 2; i++) 
	{
		for (int j=0; j < 2; j++)
		{
			output[i][j] *= _emit.emit();
		}
	}
	return;
}
