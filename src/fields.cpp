#include "fields.h"

// ==================================
// Constructors, Destructor
// ==================================
Field::Field(long _x_pts, long _y_pts, long _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag) :
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

Field::Field(SimParams &simparams) :
	x_pts(simparams.n_field_x),
	y_pts(simparams.n_field_y),
	z_pts(simparams.n_field_z),
	n_pts(simparams.n_field_x*simparams.n_field_y*simparams.n_field_z),
	x_edge_mag(simparams.radius),
	y_edge_mag(simparams.radius),
	z_edge_mag(simparams.radius)
{
	_init();
}

/* Field::Field(const Field &rhs) */
/* { */
/* 	(*this)._copy(rhs); */
/* } */

Field::~Field()
{
	delete[] x_data;
	delete[] y_data;
}

// ==================================
// Private methods
// ==================================
bool Field::_samedim(const Field &rhs)
{
	if ( ((*this).x_pts == rhs.x_pts) && ((*this).y_pts == rhs.y_pts) && ((*this).z_pts == rhs.z_pts) )
	{
		return true;
	} else {
		return false;
	}
}

/* int Field::_copy(const Field &rhs) */
/* { */
/* 	(*this).~Field(); */

/* 	(*this)._init(); */

/* 	(*this).x_data = rhs.x_data; */
/* 	(*this).y_data = rhs.y_data; */

/* 	return 0; */
/* } */

int Field::_init()
{
	x_data  = new double[n_pts];
	y_data  = new double[n_pts];
	return 0;
}

long Field::_index(long i, long j, long k)
{
	long index = i + j*x_pts + k*y_pts*x_pts;
	return index;
}

// ==================================
// Public methods
// ==================================
double &Field::x(long i, long j, long k)
{
	return x_data[_index(i, j, k)];
}

double &Field::y(long i, long j, long k)
{
	return y_data[_index(i, j, k)];
}

// ==================================
// Operators
// ==================================
Field &Field::operator+=(const Field &rhs)
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

Field &Field::operator-=(const Field &rhs)
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

const Field Field::operator+(const Field &rhs)
{
	return Field(*this) += rhs;
}

const Field Field::operator-(const Field &rhs)
{
	return Field(*this) -= rhs;
}
