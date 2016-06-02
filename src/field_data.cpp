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
	Ex_data = new double[n_pts];
	Ey_data = new double[n_pts];
	Ez_data = new double[n_pts];

	Bx_data = new double[n_pts];
	By_data = new double[n_pts];
	Bz_data = new double[n_pts];

	x_grid = new double[x_pts];
	y_grid = new double[y_pts];
	z_grid = new double[z_pts];

	for (long long i=0; i < n_pts; i++)
	{
		Ex_data[i] = double(0);
		Ey_data[i] = double(0);
		Ez_data[i] = double(0);

		Bx_data[i] = double(0);
		By_data[i] = double(0);
		Bz_data[i] = double(0);
	}

	dxdi = x_edge_mag * 2 / (x_pts-1);
	dydj = y_edge_mag * 2 / (y_pts-1);
	dzdk = z_edge_mag * 2 / (z_pts-1);

	mid_i = (x_pts-1) / 2;
	mid_j = (y_pts-1) / 2;
	mid_k = (z_pts-1) / 2;

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

	if (z_pts == 1)
	{
		z_grid[0] = 0;
	}

	splines_valid = _init_splines();

	return 0;
}

long long Field_Data::_index(int i, int j, int k) const
{
	long long index = i + x_pts*(j + y_pts*k);
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
		Ex_data[i] = rhs.Ex_data[i];
		Ey_data[i] = rhs.Ey_data[i];
		Ez_data[i] = rhs.Ez_data[i];

		Bx_data[i] = rhs.Bx_data[i];
		By_data[i] = rhs.By_data[i];
		Bz_data[i] = rhs.Bz_data[i];
	}
}

Field_Data::Field_Data(const SimParams &simparams) :
	x_pts(simparams.n_field_x),
	y_pts(simparams.n_field_y),
	z_pts(simparams.n_field_z),
	n_pts(simparams.n_field_x*simparams.n_field_y*simparams.n_field_z),
	x_edge_mag(simparams.field_trans_wind),
	y_edge_mag(simparams.field_trans_wind),
	z_edge_mag(simparams.z_end)
{
	_init();
}

Field_Data::Field_Data(int _x_pts, int _y_pts, int _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag) :
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
	delete[] Ex_data;
	delete[] Ey_data;
	delete[] Ez_data;

	delete[] Bx_data;
	delete[] By_data;
	delete[] Bz_data;

	delete[] x_grid;
	delete[] y_grid;
	delete[] z_grid;
}

// ==================================
// Public Methods:
// Facing methods
// ==================================
double &Field_Data::Ex_ind(int i, int j, int k)
{
	splines_valid = false;
	return Ex_data[_index(i, j, k)];
}

double &Field_Data::Ex_ind(long long ind)
{
	splines_valid = false;
	return Ex_data[ind];
}

double &Field_Data::Ey_ind(int i, int j, int k)
{
	splines_valid = false;
	return Ey_data[_index(i, j, k)];
}

double &Field_Data::Ey_ind(long long ind)
{
	splines_valid = false;
	return Ey_data[ind];
}

double &Field_Data::Ez_ind(int i, int j, int k)
{
	splines_valid = false;
	return Ez_data[_index(i, j, k)];
}

double &Field_Data::Ez_ind(long long ind)
{
	splines_valid = false;
	return Ez_data[ind];
}

double &Field_Data::Bx_ind(int i, int j, int k)
{
	splines_valid = false;
	return Bx_data[_index(i, j, k)];
}

double &Field_Data::Bx_ind(long long ind)
{
	splines_valid = false;
	return Bx_data[ind];
}

double &Field_Data::By_ind(int i, int j, int k)
{
	splines_valid = false;
	return By_data[_index(i, j, k)];
}

double &Field_Data::By_ind(long long ind)
{
	splines_valid = false;
	return By_data[ind];
}

double &Field_Data::Bz_ind(int i, int j, int k)
{
	splines_valid = false;
	return Bz_data[_index(i, j, k)];
}

double &Field_Data::Bz_ind(long long ind)
{
	splines_valid = false;
	return Bz_data[ind];
}

// ==================================
// Operators
// ==================================
Field_Data &Field_Data::operator+=(const Field_Data &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this).Ex_data += *rhs.Ex_data;
		*(*this).Ey_data += *rhs.Ey_data;

		*(*this).Bx_data += *rhs.Bx_data;
		*(*this).By_data += *rhs.By_data;
	} else {
		throw "Cannot add fields of different sizes";
	}
	return *this;
}

Field_Data &Field_Data::operator-=(const Field_Data &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this).Ex_data -= *rhs.Ex_data;
		*(*this).Ey_data -= *rhs.Ey_data;

		*(*this).Bx_data -= *rhs.Bx_data;
		*(*this).By_data -= *rhs.By_data;
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

