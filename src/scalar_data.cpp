#include "scalar_data.h"
#include <vector>
#include <math.h>

bool ScalarData::_samedim(const ScalarData &rhs) const
{
	if ( ((*this).x_pts == rhs.x_pts) && ((*this).y_pts == rhs.y_pts) && ((*this).z_pts == rhs.z_pts) )
	{
		return true;
	} else {
		return false;
	}
}

long long ScalarData::n_pts() const
{
	return x_pts*y_pts*z_pts;
}

int ScalarData::_init()
{
	// ==================================
	// Initialize Scalar Field
	// ==================================
	long long num_pts;
	num_pts = n_pts();
	data.resize(num_pts);

	// ==================================
	// Initialize grid indices
	// ==================================
	x_grid.resize(x_pts);
	y_grid.resize(y_pts);
	z_grid.resize(z_pts);

	// ==================================
	// Explicitly set fields to zero
	// ==================================
	for (long long i=0; i < num_pts; i++)
	{
		data[i] = double(0);
	}

	// ==================================
	// Get differential lengths
	// ==================================
	dxdi = x_edge_mag * 2 / (x_pts-1);
	dydj = y_edge_mag * 2 / (y_pts-1);
	dzdk = z_edge_mag * 2 / (z_pts-1);

	// ==================================
	// Get midpoints
	// ==================================
	mid_i = (x_pts-1) / 2;
	mid_j = (y_pts-1) / 2;
	mid_k = (z_pts-1) / 2;

	// ==================================
	// Set grid indices
	// ==================================
	for (int i=0; i < x_pts; i++)
	{
		x_grid[i] = (i-mid_i) * dxdi;
	}

	for (int j=0; j < y_pts; j++)
	{
		y_grid[j] = (j-mid_j) * dydj;
	}

	for (int k=0; k < z_pts; k++)
	{
		z_grid[k] = (k-mid_k) * dzdk;
	}

	// ==================================
	// In this case, z_grid must be set
	// explicitly
	// ==================================
	if (z_pts == 1)
	{
		z_grid[0] = 0;
	}

	return 0;
}

ScalarData::ScalarData(const ScalarData &rhs) : 
	x_pts(rhs.x_pts),
	y_pts(rhs.y_pts),
	z_pts(rhs.z_pts),
	x_edge_mag(rhs.x_edge_mag),
	y_edge_mag(rhs.y_edge_mag),
	z_edge_mag(rhs.z_edge_mag)
{
	_init();
	data = rhs.data;
}

ScalarData::ScalarData(int _x_pts, int _y_pts, int _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag) :
	x_pts(_x_pts),
	y_pts(_y_pts),
	z_pts(_z_pts),
	x_edge_mag(_x_edge_mag),
	y_edge_mag(_y_edge_mag),
	z_edge_mag(_z_edge_mag)
{
	_init();
}

ScalarData::ScalarData(const SimParams &simparams) :
	x_pts(simparams.n_field_x),
	y_pts(simparams.n_field_y),
	z_pts(simparams.n_field_z),
	x_edge_mag(simparams.field_trans_wind),
	y_edge_mag(simparams.field_trans_wind),
	z_edge_mag(simparams.z_end)
{
	_init();
}

long long ScalarData::_index(int i, int j, int k) const
{
	long long index = i + x_pts*(j + y_pts*k);
	return index;
}

long double &ScalarData::ind(int i, int j, int k)
{
	return data[_index(i, j, k)];
}

long double &ScalarData::ind(long long ind)
{
	return data[ind];
}


llong ScalarData::lt_x_ind(ldouble x)
{
	return floor(x/dxdi + mid_i);
}

llong ScalarData::lt_y_ind(ldouble y)
{
	return floor(y/dydj + mid_j);
}

llong ScalarData::lt_z_ind(ldouble z)
{
	return floor(z/dzdk + mid_k);
}

// ==================================
// Operators
// ==================================
ScalarData &ScalarData::operator=(const ScalarData &rhs)
{
	if (this != &rhs)
	{
		// Deallocate, allocate new space, copy values...
		if (
			(this->x_pts == rhs.x_pts) &
			(this->y_pts == rhs.y_pts) &
			(this->z_pts == rhs.z_pts) &
			(this->x_edge_mag == rhs.x_edge_mag) &
			(this->y_edge_mag == rhs.y_edge_mag) &
			(this->z_edge_mag == rhs.z_edge_mag)
		   )
		{
			this->data = rhs.data;
		} else {
			throw "Cannot copy fields of different sizes";
		}
	}
    	return *this;
}

ScalarData ScalarData::operator-() const
{
	ScalarData temp = *this;
	for (int i; i < n_pts(); i++)
	{
		temp.data[i] *= -1;
	}

	return temp;
}

ScalarData &ScalarData::operator+=(const ScalarData &rhs)
{
	if ( this->_samedim(rhs) )
	{
		for (int i; i < n_pts(); i++)
		{
			this->data[i] += rhs.data[i];
		}
	} else {
		throw "Cannot add fields of different sizes";
	}
	return *this;
}

ScalarData &ScalarData::operator-=(const ScalarData &rhs)
{
	if ( this->_samedim(rhs) )
	{
		for (int i; i < n_pts(); i++)
		{
			this->data[i] -= rhs.data[i];
		}
	} else {
		throw "Cannot subtract fields of different sizes";
	}
	return *this;
}

const ScalarData ScalarData::operator+(const ScalarData &rhs)
{
	return ScalarData(*this) += rhs;
}

const ScalarData ScalarData::operator-(const ScalarData &rhs)
{
	return ScalarData(*this) -= rhs;
}
