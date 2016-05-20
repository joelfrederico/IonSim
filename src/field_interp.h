#ifndef __FIELD_INTERP_H_INCLUDED__
#define __FIELD_INTERP_H_INCLUDED__

#include "field_data.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

class Field_Interp
{
	private:
		Field_Data _field;
		const gsl_interp2d_type _interptype;

		bool splines_valid;
		int z_step_current;

		double *x_data;
		double *y_data;
		double *z_data;

		gsl_interp_accel *Ex_xacc;
		gsl_interp_accel *Ex_yacc;
		gsl_interp_accel *Ey_xacc;
		gsl_interp_accel *Ey_yacc;
		gsl_interp_accel *Ez_xacc;
		gsl_interp_accel *Ez_yacc;

		/*
		gsl_spline2d *splinex;
		gsl_spline2d *spliney;
		gsl_spline2d *splinez;
		*/

		gsl_interp2d *interp2d_x;
		gsl_interp2d *interp2d_y;
		gsl_interp2d *interp2d_z;

		int _init();
		int _init_splines(int z_step);
		double field_interp(double x, double y, int z_step, int dim);
	public:
		Field_Interp(const Field_Data &field, const gsl_interp2d_type interptype);
		~Field_Interp();

		double Ex(double x, double y, int z_ind);
		double Ey(double x, double y, int z_ind);
		double Ez(double x, double y, int z_ind);
};

#endif
