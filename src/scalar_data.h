#ifndef __SCALAR_DATA_H_INCLUDED__
#define __SCALAR_DATA_H_INCLUDED__

#include "consts.h"
#include "simparams.h"

template<typename Tclass>
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
		ScalarData(const unsigned int _x_pts, const unsigned int _y_pts, const unsigned int _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag);

		// ==================================
		// Member data
		// ==================================
		long double dxdi;
		long double dydj;
		long double dzdk;

		const unsigned int x_pts;
		const unsigned int y_pts;
		const unsigned int z_pts;

		const long double x_edge_mag;
		const long double y_edge_mag;
		const long double z_edge_mag;

		std::vector<Tclass> data;
		std::vector<long double> x_grid, y_grid, z_grid;

		// ==================================
		// Calculated data
		// ==================================
		long long n_pts() const;

		// ==================================
		// Methods
		// ==================================
		unsigned long long _index(const unsigned int i, const unsigned int j, const unsigned int k) const;

		/* Tclass &ind(int i, int j, int k); */
		Tclass &ind(const unsigned int i, const unsigned int j, const unsigned int k);
		Tclass &ind(unsigned long long ind);

		int lt_x_ind_e(const ldouble x, int &ind) const;
		int lt_y_ind_e(const ldouble y, int &ind) const;
		int lt_z_ind_e(const ldouble z, int &ind) const;

		// ==================================
		// Operators
		// ==================================
		ScalarData<Tclass> &operator=(const ScalarData<Tclass> &rhs);
		ScalarData<Tclass> operator-() const;

		ScalarData<Tclass> &operator+=(const ScalarData<Tclass> &rhs);
		ScalarData<Tclass> &operator-=(const ScalarData<Tclass> &rhs);

		template <class T>
		ScalarData<Tclass> &operator=(const T rhs)
		{
			for (int i=0; i < n_pts(); i++)
			{
				this->data[i] = rhs;
			}
			return *this;
		}

		template <class T>
		ScalarData<Tclass> &operator*=(const T rhs)
		{
			for (int i=0; i < n_pts(); i++)
			{
				this->data[i] *= rhs;
			}
			return *this;
		}

		const ScalarData<Tclass> operator+(const ScalarData<Tclass> &rhs);
		const ScalarData<Tclass> operator-(const ScalarData<Tclass> &rhs);

		template <class T>
		const ScalarData<Tclass> operator*(const T rhs)
		{
			return ScalarData<Tclass>(*this) *= rhs;
		}
};

#endif
