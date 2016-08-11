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

	protected:

	public:
		double dxdi;
		double dydj;
		double dzdk;

		// ==================================
		// Constructors, Destructor
		// ==================================
		Field_Data(const Field_Data &rhs);
		Field_Data(const SimParams &simparams);
		Field_Data(int _x_pts, int _y_pts, int _z_pts, double _x_edge_mag, double _y_edge_mag, double _z_edge_mag);
		~Field_Data();

		// ==================================
		// Data members
		// ==================================
		const int x_pts;
		const int y_pts;
		const int z_pts;
		const long long n_pts;

		const double x_edge_mag;
		const double y_edge_mag;
		const double z_edge_mag;

		double *Ex_data;
		double *Ey_data;
		double *Ez_data;

		double *Bx_data;
		double *By_data;
		double *Bz_data;

		double *x_grid;
		double *y_grid;
		double *z_grid;

		// ==================================
		// Methods
		// ==================================
		long long _index(int i, int j, int k) const;

		double &Ex_ind(int i, int j, int k);
		double &Ex_ind(long long ind);

		double &Ey_ind(int i, int j, int k);
		double &Ey_ind(long long ind);

		double &Ez_ind(int i, int j, int k);
		double &Ez_ind(long long ind);

		double &Bx_ind(int i, int j, int k);
		double &Bx_ind(long long ind);

		double &By_ind(int i, int j, int k);
		double &By_ind(long long ind);

		double &Bz_ind(int i, int j, int k);
		double &Bz_ind(long long ind);

		// ==================================
		// Operators
		// ==================================
		Field_Data &operator+=(const Field_Data &rhs);
		Field_Data &operator-=(const Field_Data &rhs);

		const Field_Data operator+(const Field_Data &rhs);
		const Field_Data operator-(const Field_Data &rhs);
};

#endif
