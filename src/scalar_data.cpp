#include "scalar_data.h"
#include <vector>
#include <math.h>

template<typename T>
bool ScalarData<T>::_samedim(const ScalarData &rhs) const
{
	if ( ((*this).x_pts == rhs.x_pts) && ((*this).y_pts == rhs.y_pts) && ((*this).z_pts == rhs.z_pts) )
	{
		return true;
	} else {
		return false;
	}
}

template<typename T>
long long ScalarData<T>::n_pts() const
{
	return x_pts * y_pts * z_pts;
}

template<typename T>
int ScalarData<T>::_init()
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

template<typename T>
ScalarData<T>::ScalarData(const ScalarData &rhs) : 
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

template<typename T>
ScalarData<T>::ScalarData(const unsigned int _x_pts, const unsigned int _y_pts, const unsigned int _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag) :
	x_pts(_x_pts),
	y_pts(_y_pts),
	z_pts(_z_pts),
	x_edge_mag(_x_edge_mag),
	y_edge_mag(_y_edge_mag),
	z_edge_mag(_z_edge_mag)
{
	_init();
}

template<typename T>
ScalarData<T>::ScalarData(const SimParams &simparams) :
	x_pts(simparams.n_field_x),
	y_pts(simparams.n_field_y),
	z_pts(simparams.n_field_z),
	x_edge_mag(simparams.field_trans_wind),
	y_edge_mag(simparams.field_trans_wind),
	z_edge_mag(simparams.z_end)
{
	_init();
}

template<typename T>
unsigned long long ScalarData<T>::_index(const unsigned int i, const unsigned int j, const unsigned int k) const
{
	if ((i >= x_pts) || (j >= y_pts) || (k >= z_pts)) throw std::runtime_error("Index out of bounds!");
	unsigned long long index = i + x_pts*(j + y_pts*k);
	return index;
}

template<typename T>
T &ScalarData<T>::ind(const unsigned int i, const unsigned int j, const unsigned int k)
{
	return ind(_index(i, j, k));
}

template<typename T>
T &ScalarData<T>::ind(unsigned long long ind)
{
	if (ind > n_pts()) throw std::runtime_error("ScalarData index requested is too large.");
	return data[ind];
}

int getind_e(const ldouble x, const double delx, const double mid, const int n_pts, int &ind)
{
	ind = floor(x/delx + mid);
	if ( (0 < ind ) && (ind < n_pts-1))
	{
		return 0;
	} else {
		return -1;
	}
}

template<typename T>
int ScalarData<T>::lt_x_ind_e(const ldouble x, int &ind) const
{
	return getind_e(x, dxdi, mid_i, x_pts, ind);
	/* return floor(x/dxdi + mid_i); */
}

template<typename T>
int ScalarData<T>::lt_y_ind_e(const ldouble y, int &ind) const
{
	return getind_e(y, dydj, mid_j, y_pts, ind);
}

template<typename T>
int ScalarData<T>::lt_z_ind_e(const ldouble z, int &ind) const
{
	return getind_e(z, dzdk, 0, z_pts, ind);
}

// ==================================
// Operators
// ==================================
template<typename T>
ScalarData<T> &ScalarData<T>::operator=(const ScalarData<T> &rhs)
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
			throw std::runtime_error("Cannot copy fields of different sizes");
		}
	}
    	return *this;
}

template<typename T>
ScalarData<T> ScalarData<T>::operator-() const
{
	ScalarData temp = *this;
	for (int i=0; i < n_pts(); i++)
	{
		temp.data[i] *= -1;
	}

	return temp;
}

template<typename T>
ScalarData<T> &ScalarData<T>::operator+=(const ScalarData<T> &rhs)
{
	if ( this->_samedim(rhs) )
	{
		for (int i=0; i < n_pts(); i++)
		{
			this->data[i] += rhs.data[i];
		}
	} else {
		throw std::runtime_error("Cannot add fields of different sizes");
	}
	return *this;
}

template<typename T>
ScalarData<T> &ScalarData<T>::operator-=(const ScalarData<T> &rhs)
{
	if ( this->_samedim(rhs) )
	{
		for (int i=0; i < n_pts(); i++)
		{
			this->data[i] -= rhs.data[i];
		}
	} else {
		throw std::runtime_error("Cannot subtract fields of different sizes");
	}
	return *this;
}

template<typename T>
const ScalarData<T> ScalarData<T>::operator+(const ScalarData<T> &rhs)
{
	return ScalarData<T>(*this) += rhs;
}

template<typename T>
const ScalarData<T> ScalarData<T>::operator-(const ScalarData<T> &rhs)
{
	return ScalarData(*this) -= rhs;
}


template<typename T>
std::vector<T> ScalarData<T>::vdata() const
{
	return data;
}

template class ScalarData<long double>;
template class ScalarData<std::complex<long double>>;
