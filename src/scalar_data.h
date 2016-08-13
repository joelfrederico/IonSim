#ifndef __SCALAR_DATA_H_INCLUDED__
#define __SCALAR_DATA_H_INCLUDED__

#include "consts.h"
#include "simparams.h"

class ScalarData
{
	private:
		// ==================================
		// Private data members
		// ==================================
		double mid_i;
		double mid_j;
		double mid_k;

		// ==================================
		// Private methods
		// ==================================
		bool _samedim(const ScalarData &rhs) const;
		int _init();

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		ScalarData(const ScalarData &rhs);
		ScalarData(const SimParams &simparams);
		ScalarData(int _x_pts, int _y_pts, int _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag);

		// ==================================
		// Member data
		// ==================================
		double dxdi;
		double dydj;
		double dzdk;

		const int x_pts;
		const int y_pts;
		const int z_pts;

		const double x_edge_mag;
		const double y_edge_mag;
		const double z_edge_mag;

		ldouble_vec data;
		ldouble_vec x_grid, y_grid, z_grid;

		// ==================================
		// Calculated data
		// ==================================
		long long n_pts() const;

		// ==================================
		// Methods
		// ==================================
		long long _index(int i, int j, int k) const;

		ldouble &ind(int i, int j, int k);
		ldouble &ind(long long ind);

		int lt_x_ind_e(const ldouble x, int &ind) const;
		int lt_y_ind_e(const ldouble y, int &ind) const;
		int lt_z_ind_e(const ldouble z, int &ind) const;

		// ==================================
		// Operators
		// ==================================
		ScalarData &operator=(const ScalarData &rhs);
		ScalarData operator-() const;

		ScalarData &operator+=(const ScalarData &rhs);
		ScalarData &operator-=(const ScalarData &rhs);

		template <class T>
		ScalarData &operator=(const T rhs)
		{
			for (int i=0; i < n_pts(); i++)
			{
				this->data[i] = rhs;
			}
			return *this;
		}

		template <class T>
		ScalarData &operator*=(const T rhs)
		{
			for (int i=0; i < n_pts(); i++)
			{
				this->data[i] *= rhs;
			}
			return *this;
		}

		const ScalarData operator+(const ScalarData &rhs);
		const ScalarData operator-(const ScalarData &rhs);

		template <class T>
		const ScalarData operator*(const T rhs)
		{
			return ScalarData(*this) *= rhs;
		}
};

#endif
