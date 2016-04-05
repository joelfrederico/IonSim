#ifndef __FIELD_INTERP_H_INCLUDED__
#define __FIELD_INTERP_H_INCLUDED__

#include "field_data.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

class Field_Interp
{
	private:
		const Field_Data _field;
		const gsl_interp2d_type _interptype;

		gsl_interp_accel *Ex_xacc;
		gsl_interp_accel *Ex_yacc;
		gsl_interp_accel *Ey_xacc;
		gsl_interp_accel *Ey_yacc;
		gsl_interp_accel *Ez_xacc;
		gsl_interp_accel *Ez_yacc;

		gsl_spline2d *splinex;
		gsl_spline2d *spliney;
		gsl_spline2d *splinez;

		int _init();
	public:
		Field_Interp(const Field_Data &field, const gsl_interp2d_type interptype);
		~Field_Interp();

		double Ex(double x, double y, double z);
		double Ey(double x, double y, double z);
		double Ez(double x, double y, double z);
};

#endif
