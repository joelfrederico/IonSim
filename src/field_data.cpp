#include "field_data.h"
#include "support_func.h"
#include "simparams.h"

// ==================================
// Private Methods
// ==================================
bool Field_Data::_init_splines()
{
	return true;
}

bool Field_Data::_samedim(const Field_Data &rhs) const
{
	if ( ((*this).x_pts == rhs.x_pts) && ((*this).y_pts == rhs.y_pts) && ((*this).z_pts == rhs.z_pts) )
	{
		return true;
	} else {
		return false;
	}
}

int Field_Data::_init()
{
	x_data = new double[n_pts];
	y_data = new double[n_pts];
	z_data = new double[n_pts];

	x_grid = new double[x_pts];
	y_grid = new double[y_pts];
	z_grid = new double[z_pts];

	for (long i=0; i < n_pts; i++)
	{
		x_data[i] = double(0);
		y_data[i] = double(0);
		z_data[i] = double(0);
	}

	dxdi = x_edge_mag * 2 / (x_pts-1);
	dydj = y_edge_mag * 2 / (y_pts-1);
	dzdk = z_edge_mag * 2 / (z_pts-1);

	mid_i = (x_pts-1) / 2;
	mid_j = (y_pts-1) / 2;
	mid_k = (z_pts-1) / 2;

	for (long i=0; i < x_pts; i++)
	{
		x_grid[i] = (i-mid_i) * dxdi;
	}

	for (long j=0; j < y_pts; j++)
	{
		y_grid[j] = (j-mid_j) * dydj;
	}

	for (long k=0; k < z_pts; k++)
	{
		z_grid[k] = (k-mid_k) * dzdk;
	}

	if (z_pts == 1)
	{
		z_grid[0] = 0;
	}

	splines_valid = _init_splines();

	return 0;
}

long Field_Data::_index(long i, long j, long k) const
{
	long index = i + x_pts*(j + y_pts*k);
	return index;
}

// ==================================
// Public Methods:
// Constructors/Destructors
// ==================================
Field_Data::Field_Data(const Field_Data &rhs) : 
	x_pts(rhs.x_pts),
	y_pts(rhs.y_pts),
	z_pts(rhs.z_pts),
	n_pts(rhs.n_pts),
	x_edge_mag(rhs.x_edge_mag),
	y_edge_mag(rhs.y_edge_mag),
	z_edge_mag(rhs.z_edge_mag)
{
	_init();
	for (int i=0; i < n_pts; i++)
	{
		x_data[i] = rhs.x_data[i];
		y_data[i] = rhs.y_data[i];
		z_data[i] = rhs.z_data[i];
	}
}

Field_Data::Field_Data(const SimParams &simparams) :
	x_pts(simparams.n_field_x),
	y_pts(simparams.n_field_y),
	z_pts(simparams.n_field_z),
	n_pts(simparams.n_field_x*simparams.n_field_y*simparams.n_field_z),
	x_edge_mag(simparams.radius * 3),
	y_edge_mag(simparams.radius * 3),
	z_edge_mag(simparams.length)
{
	_init();
}

Field_Data::Field_Data(long _x_pts, long _y_pts, long _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag) :
	x_pts(_x_pts),
	y_pts(_y_pts),
	z_pts(_z_pts),
	n_pts(_x_pts*_y_pts*_z_pts),
	x_edge_mag(_x_edge_mag),
	y_edge_mag(_y_edge_mag),
	z_edge_mag(_z_edge_mag)
{
	_init();
}

Field_Data::~Field_Data()
{
	delete[] x_data;
	delete[] y_data;
	delete[] z_data;

	delete[] x_grid;
	delete[] y_grid;
	delete[] z_grid;
}

// ==================================
// Public Methods:
// Facing methods
// ==================================
double &Field_Data::Ex_ind(long i, long j, long k)
{
	splines_valid = false;
	return x_data[_index(i, j, k)];
}

double &Field_Data::Ex_ind(long ind)
{
	splines_valid = false;
	return x_data[ind];
}

/* double Field_Data::Ex_ind(long i, long j, long k) const */
/* { */
/* 	return x_data[_index(i, j, k)]; */
/* } */

double &Field_Data::Ey_ind(long i, long j, long k)
{
	splines_valid = false;
	return y_data[_index(i, j, k)];
}

double &Field_Data::Ey_ind(long ind)
{
	splines_valid = false;
	return y_data[ind];
}

double Field_Data::Ey_ind(long i, long j, long k) const
{
	return y_data[_index(i, j, k)];
}

double &Field_Data::Ez_ind(long i, long j, long k)
{
	splines_valid = false;
	return z_data[_index(i, j, k)];
}

double &Field_Data::Ez_ind(long ind)
{
	splines_valid = false;
	return z_data[ind];
}

double Field_Data::Ez_ind(long i, long j, long k) const
{
	return z_data[_index(i, j, k)];
}


// ==================================
// Operators
// ==================================
Field_Data &Field_Data::operator+=(const Field_Data &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this).x_data += *rhs.x_data;
		*(*this).y_data += *rhs.y_data;
	} else {
		throw "Cannot add fields of different sizes";
	}
	return *this;
}

Field_Data &Field_Data::operator-=(const Field_Data &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this).x_data -= *rhs.x_data;
		*(*this).y_data -= *rhs.y_data;
	} else {
		throw "Cannot subtract fields of different sizes";
	}
	return *this;
}

const Field_Data Field_Data::operator+(const Field_Data &rhs)
{
	return Field_Data(*this) += rhs;
}

const Field_Data Field_Data::operator-(const Field_Data &rhs)
{
	return Field_Data(*this) -= rhs;
}

