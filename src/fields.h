#ifndef __FIELDS_H_INCLUDED__
#define __FIELDS_H_INCLUDED__

#include "consts.h"
#include "parts.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <mpi.h>

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

		// ==================================
		// Private methods
		// ==================================
		bool _samedim(const Field &rhs);
		int _init();
		long _index(long i, long j) const;

		int _array_alloc(double ** (&out), double *data);

		/* gsl_spline2d *splinex; */
		/* gsl_spline2d *spliney; */
		/* gsl_interp_accel *Ex_xacc; */
		/* gsl_interp_accel *Ex_yacc; */
		/* gsl_interp_accel *Ey_xacc; */
		/* gsl_interp_accel *Ey_yacc; */
		bool splines_valid;

		bool _init_splines();

		double *x_grid, *y_grid;
	public:
		/* const gsl_interp2d_type *T; */
		// ==================================
		// Constructors, Destructor
		// ==================================
		Field(long _x_pts, long _y_pts, long _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag);
		Field(const SimParams &simparams);
		Field(const Field &rhs);
		~Field();

		// ==================================
		// Data members
		// ==================================
		const long x_pts;
		const long y_pts;
		const long z_pts;

		const double x_edge_mag;
		const double y_edge_mag;
		const double z_edge_mag;

		double *x_data;
		double *y_data;
		double *z_data;

		// ==================================
		// Methods
		// ==================================
		double &Ex_ind(long i, long j, long k);
		double Ex_ind(long i, long j, long k) const;
		double Ex(double x, double y, double z);

		double &Ey_ind(long i, long j, long k);
		double Ey_ind(long i, long j, long k) const;
		double Ey(double x, double y, double z);

		double &Ez_ind(long i, long j, long k);
		double Ez_ind(long i, long j, long k) const;
		double Ez(double x, double y, double z);

		double i_to_x(long i);
		double i(double _x);

		double j_to_y(long j);
		double j(double _y);

		double k_to_z(long k);
		double k(double _z);

		int dump_serial(std::string const &filename, long step);

		int recv_field_others();
		int recv_field(int sender_id);
		int send_field(int dest_id);

		int get_interp(Field &field_in);

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
				*(*this).z_data *= rhs;
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
				*(*this).z_data /= rhs;
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
