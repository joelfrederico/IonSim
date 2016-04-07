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

	x_data = new double[_field.x_pts*_field.y_pts];
	y_data = new double[_field.x_pts*_field.y_pts];
	z_data = new double[_field.x_pts*_field.y_pts];

	z_step_current = -1;

	splines_valid = false;
	return 0;
}

int Field_Interp::_init_splines(int z_step)
{
	gsl_spline2d_init(splinex, _field.x_grid, _field.y_grid, x_data, _field.x_pts, _field.y_pts);
	gsl_spline2d_init(spliney, _field.x_grid, _field.y_grid, y_data, _field.x_pts, _field.y_pts);
	gsl_spline2d_init(splinez, _field.x_grid, _field.y_grid, z_data, _field.x_pts, _field.y_pts);

	for (int i=0; i < _field.x_pts; i++)
	{
		for (int j=0; j < _field.y_pts; j++)
		{
			gsl_spline2d_set(splinex, x_data, i, j, _field.Ex_ind(i, j, z_step));
			gsl_spline2d_set(spliney, y_data, i, j, _field.Ey_ind(i, j, z_step));
			gsl_spline2d_set(splinez, z_data, i, j, _field.Ez_ind(i, j, z_step));
		}
	}

	splines_valid = true;

	return 0;
}

Field_Interp::~Field_Interp()
{
	delete x_data;
	delete y_data;
	delete z_data;

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

double Field_Interp::Ex(double x, double y, int z_step)
{
	return field_interp(x, y, z_step, 0);
}

double Field_Interp::Ey(double x, double y, int z_step)
{
	return field_interp(x, y, z_step, 1);
}

double Field_Interp::Ez(double x, double y, int z_step)
{
	return field_interp(x, y, z_step, 2);
}

double Field_Interp::field_interp(double x, double y, int z_step, int dim)
{
	int err;
	double out;

	gsl_spline2d *spline;
	gsl_interp_accel *xacc;
	gsl_interp_accel *yacc;

	if (z_step != z_step_current)
	{
		_init_splines(z_step);
	}

	switch (dim)
	{
		case 0:
			xacc = Ex_xacc;
			yacc = Ex_yacc;
			spline = splinex;
			break;
		case 1:
			xacc = Ey_xacc;
			yacc = Ey_yacc;
			spline = spliney;
			break;
		case 2:
			xacc = Ez_xacc;
			yacc = Ez_yacc;
			spline = splinez;
			break;
	}


	err = gsl_spline2d_eval_e(spline, x, y, xacc, yacc, &out);

	return out;
}
