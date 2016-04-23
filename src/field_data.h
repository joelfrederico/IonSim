#ifndef __FIELD_DATA_H_INCLUDED__
#define __FIELD_DATA_H_INCLUDED__

#include "consts.h"
#include <mpi.h>
#include "simparams.h"

class Field_Data
{
	private:
		// ==================================
		// Private data members
		// ==================================
		double dxdi;
		double dydj;

		double mid_i;
		double mid_j;
		double mid_k;

		bool splines_valid;

		// ==================================
		// Private methods
		// ==================================
		bool _init_splines();
		bool _samedim(const Field_Data &rhs) const;
		int _init();
		long _index(long i, long j, long k) const;

	protected:

	public:
		double dzdk;

		// ==================================
		// Constructors, Destructor
		// ==================================
		Field_Data(const Field_Data &rhs);
		Field_Data(const SimParams &simparams);
		Field_Data(long _x_pts, long _y_pts, long _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag);
		~Field_Data();

		// ==================================
		// Data members
		// ==================================
		const long x_pts;
		const long y_pts;
		const long z_pts;
		const long n_pts;

		const double x_edge_mag;
		const double y_edge_mag;
		const double z_edge_mag;

		double *x_data;
		double *y_data;
		double *z_data;

		double *x_grid;
		double *y_grid;
		double *z_grid;

		// ==================================
		// Methods
		// ==================================
		double &Ex_ind(long i, long j, long k);
		double &Ex_ind(long ind);
		double Ex_ind(long i, long j, long k) const;

		double &Ey_ind(long i, long j, long k);
		double &Ey_ind(long ind);
		double Ey_ind(long i, long j, long k) const;

		double &Ez_ind(long i, long j, long k);
		double &Ez_ind(long ind);
		double Ez_ind(long i, long j, long k) const;

		// ==================================
		// Operators
		// ==================================
		Field_Data &operator+=(const Field_Data &rhs);
		Field_Data &operator-=(const Field_Data &rhs);

		const Field_Data operator+(const Field_Data &rhs);
		const Field_Data operator-(const Field_Data &rhs);
};

#endif
