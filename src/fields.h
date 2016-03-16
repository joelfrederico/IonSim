#ifndef __FIELDS_H_INCLUDED__
#define __FIELDS_H_INCLUDED__

#include "consts.h"
#include "parts.h"

class Field;

class Field
{
	private:
		// ==================================
		// Data members
		// ==================================
		long _n_pts;

		// ==================================
		// Methods
		// ==================================
		bool _samedim(const Field &rhs);
		int _copy(const Field &rhs);
		int _init(long x_pts, long y_pts, long z_pts);
		long _index(long i, long j, long k);

	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		Field(long _x_pts, long _y_pts, long _z_pts);
		Field(SimParams &simparams);
		Field(const Field &rhs);
		~Field();

		// ==================================
		// Data members
		// ==================================
		long x_pts;
		long y_pts;
		long z_pts;

		double *x_data;
		double *y_data;

		// ==================================
		// Methods
		// ==================================
		double &x(long i, long j, long k);
		double &y(long i, long j, long k);

		// ==================================
		// Operators
		// ==================================
		Field &operator=(const Field &rhs);
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
