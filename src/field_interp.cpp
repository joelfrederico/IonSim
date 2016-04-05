#include "field_interp.h"
#include "field_data.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

Field_Interp::Field_Interp(const Field_Data &field, const gsl_interp2d_type interptype) : 
	_field(field),
	_interptype(interptype)
{
	_init();
}

int Field_Interp::_init()
{
	Ex_xacc = gsl_interp_accel_alloc();
	Ex_yacc = gsl_interp_accel_alloc();
	Ey_xacc = gsl_interp_accel_alloc();
	Ey_yacc = gsl_interp_accel_alloc();
	Ez_xacc = gsl_interp_accel_alloc();
	Ez_yacc = gsl_interp_accel_alloc();

	splinex = gsl_spline2d_alloc(&_interptype, _field.x_pts, _field.y_pts);
	spliney = gsl_spline2d_alloc(&_interptype, _field.x_pts, _field.y_pts);
	splinez = gsl_spline2d_alloc(&_interptype, _field.x_pts, _field.y_pts);

	gsl_spline2d_init(splinex, _field.x_grid, _field.y_grid, _field.x_data, _field.x_pts, _field.y_pts);
	gsl_spline2d_init(spliney, _field.x_grid, _field.y_grid, _field.y_data, _field.x_pts, _field.y_pts);
	gsl_spline2d_init(splinez, _field.x_grid, _field.y_grid, _field.z_data, _field.x_pts, _field.y_pts);

	return 0;
}

Field_Interp::~Field_Interp()
{
	gsl_spline2d_free(splinex);
	gsl_spline2d_free(spliney);
	gsl_spline2d_free(splinez);

	gsl_interp_accel_free(Ex_xacc);
	gsl_interp_accel_free(Ex_yacc);
	gsl_interp_accel_free(Ey_xacc);
	gsl_interp_accel_free(Ey_yacc);
	gsl_interp_accel_free(Ez_xacc);
	gsl_interp_accel_free(Ez_yacc);

}

double Field_Interp::Ex(double x, double y, double z)
{
	int err;
	double out;
	err = gsl_spline2d_eval_e(splinex, x, y, Ex_xacc, Ex_yacc, &out);
	return out;
}

double Field_Interp::Ey(double x, double y, double z)
{
	int err;
	double out;
	err = gsl_spline2d_eval_e(spliney, x, y, Ey_xacc, Ey_yacc, &out);
	return out;
}

double Field_Interp::Ez(double x, double y, double z)
{
	int err;
	double out;
	err = gsl_spline2d_eval_e(splinez, x, y, Ez_xacc, Ez_yacc, &out);
	return out;
}
