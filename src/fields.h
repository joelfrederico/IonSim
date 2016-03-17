#ifndef __FIELDS_H_INCLUDED__
#define __FIELDS_H_INCLUDED__

#include "consts.h"
#include "parts.h"

class Field;

class Field
{
	private:
		// ==================================
		// Private data members
		// ==================================
		const long n_pts;
		double dxdi;
		double dydj;
		double dzdk;
		double mid_i;
		double mid_j;
		double mid_k;

		// ==================================
		// Private methods
		// ==================================
		bool _samedim(const Field &rhs);
		/* int _copy(const Field &rhs); */
		int _init();
		long _index(long i, long j, long k);

		int _array_alloc(double ** (&out), long k, double *data);

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		Field(long _x_pts, long _y_pts, long _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_length);
		Field(SimParams &simparams);
		~Field();

		// ==================================
		// Data members
		// ==================================
		const long x_pts;
		const long y_pts;
		const long z_pts;

		const double x_edge_mag;
		const double y_edge_mag;
		const double z_length;

		double *x_data;
		double *y_data;

		// ==================================
		// Methods
		// ==================================
		double &x(long i, long j, long k);
		double &y(long i, long j, long k);
		double i_to_x(long i);
		double j_to_y(long j);
		double k_to_z(long k);
		double i(double _x, double _y);
		double j(double _x, double _y);

		int x_array_alloc(double ** (&out), long k);
		int y_array_alloc(double ** (&out), long k);

		// ==================================
		// Operators
		// ==================================
		Field &operator+=(const Field &rhs);
		Field &operator-=(const Field &rhs);
		template <class T>
		Field &operator*=(const T rhs)
		{
			if ( (*this)._samedim(rhs) )
			{
				*(*this).x_data *= rhs;
				*(*this).y_data *= rhs;
			} else {
				throw "Cannot subtract fields of different sizes";
			}
			return *this;
		}
		template <class T>
		Field &operator/=(const T rhs)
		{
			if ( (*this)._samedim(rhs) )
			{
				*(*this).x_data /= rhs;
				*(*this).y_data /= rhs;
			} else {
				throw "Cannot subtract fields of different sizes";
			}
			return *this;
		}

		const Field operator+(const Field &rhs);
		const Field operator-(const Field &rhs);

		template <class T>
		const Field operator*(T rhs)
		{
			return Field(*this) *= rhs;
		}
		template <class T>
		const Field operator/(T rhs)
		{
			return Field(*this) /= rhs;
		}
};

#endif
