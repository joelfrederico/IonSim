#include "fields.h"

Field::Field(long x_pts, long y_pts)
{
	(*this)._init(x_pts, y_pts);
}

Field::~Field()
{
	delete[] _x_data;
	delete[] _y_data;
}

int Field::_init(long x_pts, long y_pts)
{
	_x_pts = x_pts;
	_y_pts = y_pts;
	_n_pts = x_pts * y_pts;
	_x_data  = new double[_n_pts];
	_y_data  = new double[_n_pts];
	return 0;
}

Field::Field(const Field &rhs)
{
	(*this)._copy(rhs);
}

int Field::_copy(const Field &rhs)
{
	(*this).~Field();

	(*this)._init(rhs._x_pts, rhs._y_pts);

	(*this)._x_data = rhs._x_data;
	(*this)._y_data = rhs._y_data;

	return 0;
}

bool Field::_samedim(const Field &rhs)
{
	if ( ((*this)._x_pts == rhs._x_pts) && ((*this)._y_pts == rhs._y_pts) )
	{
		return true;
	} else {
		return false;
	}
}

Field &Field::operator=(const Field &rhs)
{
	if (this != &rhs)
	{
		(*this)._copy(rhs);
	}
	return *this;
}

Field &Field::operator+=(const Field &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this)._x_data += *rhs._x_data;
		*(*this)._y_data += *rhs._y_data;
	} else {
		throw "Cannot add fields of different sizes";
	}
	return *this;
}

Field &Field::operator-=(const Field &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this)._x_data -= *rhs._x_data;
		*(*this)._y_data -= *rhs._y_data;
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
