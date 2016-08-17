#include "scalar_data.h"
#include <vector>
#include <math.h>

// ==================================
// Private
// ==================================
template<typename T>
bool ScalarData<T>::_samedim(const ScalarData<T> &rhs) const
{
	bool out = true;

	// Must have same rank
	if (_x_pts.size() != rhs._x_pts.size()) return false;

	// Must have same dims
	for (typename decltype(_x_pts)::size_type i=0; i<_x_pts.size(); i++)
	{
		if (_x_pts[i] != rhs._x_pts[i])
		{
			return false;
		}
	}
	return true;
}

template<typename T>
int ScalarData<T>::_init(
		const std::vector<typename std::vector<T>::size_type> x_pts,
		const std::vector<long double> edge_mag
	)
{
	// ==================================
	// Check if rank makes sense
	// ==================================
	if (x_pts.size() != edge_mag.size()) throw std::runtime_error("Inconsistent rank!");

	typename decltype(data)::size_type _n_pts = 1;

	const typename decltype(x_pts)::size_type rank = x_pts.size();

	// ==================================
	// Save private members first
	// ==================================
	// Initialization parameters
	_x_pts    = x_pts;
	_edge_mag = edge_mag;

	// Metadata parameters
	mid.resize(rank);
	_dxdi.resize(rank);
	_grid.resize(rank);
	for (typename decltype(x_pts)::size_type i=0; i<rank; i++)
	{
		if (x_pts[i] <= 1) throw std::runtime_error("Makes no sense having a dimension with no extent (x_pts=1).");
		// Obtain product
		_n_pts *= x_pts[i];
		// Calculate midpoint
		mid[i] = (x_pts[i]-1) / 2;
		// Calculate differential length
		_dxdi[i] = edge_mag[i] * 2 / (x_pts[i]-1);
		// Calculate grid vector
		_grid[i].resize(x_pts[i]);
		for (typename decltype(x_pts)::value_type j=0; j<x_pts[i]; j++)
		{
			_grid[i][j] = (j-mid[i]) * _dxdi[i];
		}
	}

	// Data storage
	data = decltype(data)(_n_pts, 0);

	return 0;
}

// ==================================
// Public
// ==================================

template<typename T>
ScalarData<T>::ScalarData()
{
	decltype(_x_pts) x_pts;
	decltype(_edge_mag) edge_mag;
	_init(x_pts, edge_mag);
}

template<typename T>
ScalarData<T>::ScalarData(const SimParams &simparams)
{
	decltype(_x_pts) x_pts = {
		simparams.n_field_x,
		simparams.n_field_y,
		simparams.n_field_z
		};
	decltype(_edge_mag) edge_mag = {
		simparams.field_trans_wind,
		simparams.field_trans_wind,
		simparams.z_end
		};
	_init(x_pts, edge_mag);
}

template<typename T>
ScalarData<T>::ScalarData(const decltype(_x_pts) x_pts, const decltype(_edge_mag) edge_mag)
{
	_init(x_pts, edge_mag);
}

template<typename T>
auto ScalarData<T>::dxdi(const decltype(_dxdi)::size_type i) const -> typename decltype(_dxdi)::value_type
{
	return _dxdi[i];
}

template<typename T>
auto ScalarData<T>::x_pts(const decltype(_x_pts)::size_type i) const -> typename decltype(_x_pts)::value_type
{
	return _x_pts[i];
}

template<typename T>
auto ScalarData<T>::n_pts() const -> typename decltype(data)::size_type
{
	return data.size();
}

template<typename T>
auto ScalarData<T>::edge_mag(const decltype(_edge_mag)::size_type i) const -> typename decltype(_edge_mag)::value_type
{
	return _edge_mag[i];
}

template<typename T>
auto ScalarData<T>::grid(const decltype(_grid)::size_type i) const -> typename decltype(_grid)::value_type
{
	return _grid[i];
}

template<typename T>
auto ScalarData<T>::vdata() const -> decltype(data)
{
	return data;
}

template<typename T>
T &ScalarData<T>::ind(typename decltype(data)::size_type i)
{
	if (i >= data.size()) throw std::runtime_error("Index is too large for data vector.");
	return data[i];
}

template<typename T>
int ScalarData<T>::lt_x_ind_e(const typename decltype(_x_pts)::size_type i, const long double x, decltype(_x_pts)::value_type &ind) const
{
	ind = floor(x/dxdi(i) + mid[i]);
	if ( (0 < ind ) && (ind < n_pts()-1))
	{
		return 0;
	} else {
		return -1;
	}
}

// ==================================
// Operators
// ==================================
/* template<typename T> */
/* ScalarData<T> &ScalarData<T>::operator=(const ScalarData<T> &rhs) */
/* { */
/* 	if (this != &rhs) */
/* 	{ */
/* 		for (typename decltype(_x_pts)::size_type i=0; i<_x_pts.size(); i++) */
/* 		{ */
/* 			if (this->x_pts(i) != rhs.x_pts(i)) throw std::runtime_error("Cannot copy, fields different sizes."); */
/* 			if (this->edge_mag(i) != rhs.edge_mag(i)) throw std::runtime_error("Cannot copy, field grids different."); */
/* 		} */
/* 		// Deallocate, allocate new space, copy values... */
/* 		this->data = rhs.data; */
/* 	} */
/*     	return *this; */
/* } */

template<typename T>
ScalarData<T> &ScalarData<T>::operator=(const T rhs)
{
	for (int i=0; i < n_pts(); i++)
	{
		this->data[i] = rhs;
	}
	return *this;
}

template<typename T>
ScalarData<T> &ScalarData<T>::operator*=(const T rhs)
{
	for (typename decltype(data)::size_type i=0; i < n_pts(); i++)
	{
		this->data[i] *= rhs;
	}
	return *this;
}

template<typename T>
ScalarData<T> ScalarData<T>::operator-() const
{
	ScalarData temp = *this;
	for (typename decltype(data)::size_type i=0; i < n_pts(); i++)
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
		for (typename decltype(data)::size_type i=0; i < n_pts(); i++)
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
		for (typename decltype(data)::size_type i=0; i < n_pts(); i++)
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

template class ScalarData<long double>;
template class ScalarData<std::complex<long double>>;
