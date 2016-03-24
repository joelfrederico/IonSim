#include "fields.h"
#include "support_func.h"
#include <gsl/gsl_interp2d.h>

// ==================================
// Constructors, Destructor
// ==================================
Field::Field(long _x_pts, long _y_pts, double _x_edge_mag, double _y_edge_mag) :
	x_pts(_x_pts),
	y_pts(_y_pts),
	n_pts(_x_pts*_y_pts),
	x_edge_mag(_x_edge_mag),
	y_edge_mag(_y_edge_mag),
	T(gsl_interp2d_bicubic)
{
	_init();
}

Field::Field(SimParams &simparams) :
	x_pts(simparams.n_field_x),
	y_pts(simparams.n_field_y),
	n_pts(simparams.n_field_x*simparams.n_field_y),
	x_edge_mag(simparams.radius),
	y_edge_mag(simparams.radius),
	T(gsl_interp2d_bicubic)
{
	_init();
}

Field::~Field()
{
	delete[] x_data;
	delete[] y_data;
	delete[] x_grid;
	delete[] y_grid;
	gsl_spline2d_free(splinex);
	gsl_spline2d_free(spliney);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc);
}

// ==================================
// Private methods
// ==================================
bool Field::_samedim(const Field &rhs)
{
	if ( ((*this).x_pts == rhs.x_pts) && ((*this).y_pts == rhs.y_pts) )
	{
		return true;
	} else {
		return false;
	}
}

int Field::_init()
{
	x_data = new double[n_pts];
	y_data = new double[n_pts];
	x_grid = new double[x_pts];
	y_grid = new double[y_pts];

	for (long i=0; i < n_pts; i++)
	{
		x_data[i] = double(0);
		y_data[i] = double(0);
	}

	splinex = gsl_spline2d_alloc(T, x_pts, y_pts);
	spliney = gsl_spline2d_alloc(T, x_pts, y_pts);
	xacc    = gsl_interp_accel_alloc();
	yacc    = gsl_interp_accel_alloc();

	dxdi = x_edge_mag * 2 / (x_pts-1);
	dydj = y_edge_mag * 2 / (y_pts-1);

	mid_i = x_pts / 2;
	mid_j = y_pts / 2;
	for (long i=0; i < x_pts; i++)
	{
		x_grid[i] = i_to_x(i);
		/* std::cout << x_grid[i] << std::endl; */
	}
	for (long j=0; j < y_pts; j++)
	{
		y_grid[j] = j_to_y(j);
		/* std::cout << y_grid[j] << std::endl; */
	}

	splines_valid = _init_splines();

	return 0;
}

bool Field::_init_splines()
{
	gsl_spline2d_init(splinex, x_grid, y_grid, x_data, x_pts, y_pts);
	gsl_spline2d_init(spliney, x_grid, y_grid, y_data, x_pts, y_pts);
	
	return true;
}

long Field::_index(long i, long j) const
{
	long index = i + j*x_pts;
	return index;
}

double Field::i_to_x(long i)
{
	return (i-mid_i) * dxdi;
}

double Field::j_to_y(long j)
{
	return (j-mid_j) * dydj;
}

int Field::_array_alloc(double ** (&out), double* data)
{
	out = ionsim::alloc_2d_array(x_pts, y_pts);
	for (int i=0; i < x_pts; i++)
	{
		for (int j=0; j < y_pts; j++)
		{
			out[i][j] = data[_index(i, j)];
		}
	}
	return x_pts;
}

int Field::x_array_alloc(double ** (&out), long k)
{
	return _array_alloc(out, x_data);
}

int Field::y_array_alloc(double ** (&out), long k)
{
	return _array_alloc(out, y_data);
}

// ==================================
// Public methods
// ==================================
double &Field::Ex_ind(long i, long j)
{
	splines_valid = false;
	return x_data[_index(i, j)];
}

double Field::Ex_ind(long i, long j) const
{
	return x_data[_index(i, j)];
}

double &Field::Ey_ind(long i, long j)
{
	splines_valid = false;
	return y_data[_index(i, j)];
}

double Field::Ey_ind(long i, long j) const
{
	return y_data[_index(i, j)];
}

double Field::Ex(double x, double y)
{
	if (!splines_valid)
	{
		splines_valid = _init_splines();
	}
	return gsl_spline2d_eval(splinex, x, y, xacc, yacc);
}

double Field::Ey(double x, double y)
{
	if (!splines_valid)
	{
		splines_valid = _init_splines();
	}
	return gsl_spline2d_eval(spliney, x, y, xacc, yacc);
}

double Field::i(double _x, double _y)
{
	return _x / dxdi + mid_i;
}

double Field::j(double _x, double _y)
{
	return _y / dydj + mid_j;
}

int Field::dump_serial(std::string const &filename, long step)
{
	/* ionsim::serial_dump(filename, step, std::string const &group, std::string const &dataset, MPI::Intracomm &comm, const Field &field) */
	std::string group;
	group = "field";

	ionsim::dump(filename, step, group, *this);
	return 0;
}

// ==================================
// Operators
// ==================================
Field &Field::operator+=(const Field &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this).x_data += *rhs.x_data;
		*(*this).y_data += *rhs.y_data;
	} else {
		throw "Cannot add fields of different sizes";
	}
	return *this;
}

Field &Field::operator-=(const Field &rhs)
{
	if ( (*this)._samedim(rhs) )
	{
		*(*this).x_data -= *rhs.x_data;
		*(*this).y_data -= *rhs.y_data;
	} else {
		throw "Cannot subtract fields of different sizes";
	}
	return *this;
}

const Field Field::operator+(const Field &rhs)
{
	return Field(*this) += rhs;
}

const Field Field::operator-(const Field &rhs)
{
	return Field(*this) -= rhs;
}
