#include "fields.h"
#include "support_func.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>

const gsl_interp2d_type *INTERPTYPE = gsl_interp2d_bicubic;

// ==================================
// Constructors, Destructor
// ==================================
Field::Field(long _x_pts, long _y_pts, double _x_edge_mag, double _y_edge_mag) :
	x_pts(_x_pts),
	y_pts(_y_pts),
	n_pts(_x_pts*_y_pts),
	x_edge_mag(_x_edge_mag),
	y_edge_mag(_y_edge_mag),
	T(INTERPTYPE)
{
	_init();
}

Field::Field(const SimParams &simparams) :
	x_pts(simparams.n_field_x),
	y_pts(simparams.n_field_y),
	n_pts(simparams.n_field_x*simparams.n_field_y),
	x_edge_mag(simparams.radius * 3),
	y_edge_mag(simparams.radius * 3),
	T(INTERPTYPE)
{
	_init();
}

Field::~Field()
{
	delete[] x_data;
	delete[] y_data;
	delete[] x_grid;
	delete[] y_grid;
	if (splines_valid)
	{
		gsl_spline2d_free(splinex);
		gsl_spline2d_free(spliney);
	}
	gsl_interp_accel_free(Ex_xacc);
	gsl_interp_accel_free(Ex_yacc);
	gsl_interp_accel_free(Ey_xacc);
	gsl_interp_accel_free(Ey_yacc);
}

Field::Field(const Field &rhs) : 
	x_pts(rhs.x_pts),
	y_pts(rhs.y_pts),
	n_pts(rhs.n_pts),
	x_edge_mag(rhs.x_edge_mag),
	y_edge_mag(rhs.y_edge_mag),
	T(rhs.T)
{
	_init();
	for (int i=0; i < n_pts; i++)
	{
		x_data[i] = rhs.x_data[i];
		y_data[i] = rhs.y_data[i];
	}
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

	Ex_xacc    = gsl_interp_accel_alloc();
	Ex_yacc    = gsl_interp_accel_alloc();
	Ey_xacc    = gsl_interp_accel_alloc();
	Ey_yacc    = gsl_interp_accel_alloc();

	dxdi = x_edge_mag * 2 / (x_pts-1);
	dydj = y_edge_mag * 2 / (y_pts-1);

	mid_i = (x_pts-1) / 2;
	mid_j = (y_pts-1) / 2;
	for (long i=0; i < x_pts; i++)
	{
		x_grid[i] = i_to_x(i);
		/* std::cout << "x[" << i << "]=" << x_grid[i] << std::endl; */
	}
	for (long j=0; j < y_pts; j++)
	{
		y_grid[j] = j_to_y(j);
		/* std::cout << "y[" << j << "]=" << y_grid[j] << std::endl; */
		/* std::cout << y_grid[j] << std::endl; */
	}

	splines_valid = _init_splines();

	return 0;
}

bool Field::_init_splines()
{
	int status;
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
	int err;
	double out;
	if (!splines_valid)
	{
		splines_valid = _init_splines();
		if (!splines_valid) printf("Something still wrong...");
	}
	err = gsl_spline2d_eval_e(splinex, x, y, Ex_xacc, Ex_yacc, &out);
	if (err != 0) printf("Eval error: %d\n", err);
	return out;
}

double Field::Ey(double x, double y)
{
	if (!splines_valid)
	{
		splines_valid = _init_splines();
	}

	return gsl_spline2d_eval(spliney, x, y, Ey_xacc, Ey_yacc);
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
	/* ionsim::dump_serial(filename, step, group, *this); */
	// ==================================
	// Initialize all variables
	// ==================================
	double *xbuf;
	double *ybuf;
	std::string x_dataset, y_dataset;
	long x_len = x_pts;
	long y_len = y_pts;

	int n_write = ionsim::MAX_N_WRITE;

	hid_t file_id;
	hid_t step_group_id, group_id;
	hid_t x_dataspace_id, y_dataspace_id;
	hid_t x_dataset_id, y_dataset_id;
	hid_t x_memspace_id, y_memspace_id;
	hsize_t count[2];
	hsize_t offset[2];
	std::string temp_str;
	herr_t status;
	
	// ==================================
	// Open/create file
	// ==================================
	file_id = ionsim::open_file(filename);
	
	// ==================================
	// Access or create a new group
	// ==================================
	step_group_id    = ionsim::group_step_access(file_id, step);

	group_id         = ionsim::group_access(step_group_id, "field");

	// ==================================
	// Create dataset
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	count[0] = x_len;
	count[1] = y_len;
	x_dataset = "Ex";
	x_dataset_id = ionsim::dataset_create(group_id, x_dataspace_id, count, x_dataset); // Updates dataspace_id
	x_memspace_id = H5Screate_simple(2, count, NULL);

	count[0] = x_len;
	count[1] = y_len;
	y_dataset = "Ey";
	y_dataset_id = ionsim::dataset_create(group_id, y_dataspace_id, count, y_dataset); // Updates dataspace_id
	y_memspace_id = H5Screate_simple(2, count, NULL);

	// ==================================
	// Write dataset rows
	// ==================================
	xbuf = new double[x_len*y_len];
	ybuf = new double[x_len*y_len];
	for (int i=0; i < x_len; i++)
	{
		for (int j=0; j < y_len; j++)
		{
			xbuf[i + j*x_len] = Ex_ind(i, j);
			ybuf[i + j*x_len] = Ey_ind(i, j);
		}
	}
	H5Dwrite(x_dataset_id, H5T_NATIVE_DOUBLE, x_memspace_id, x_dataspace_id, H5P_DEFAULT, xbuf);
	H5Dwrite(y_dataset_id, H5T_NATIVE_DOUBLE, y_memspace_id, y_dataspace_id, H5P_DEFAULT, ybuf);

	delete [] xbuf;
	delete [] ybuf;

	ionsim::writeattribute(group_id , "n_pts"      , n_pts);
	ionsim::writeattribute(group_id , "dxdi"       , dxdi);
	ionsim::writeattribute(group_id , "dydj"       , dydj);
	ionsim::writeattribute(group_id , "mid_i"      , mid_i);
	ionsim::writeattribute(group_id , "mid_j"      , mid_j);
	ionsim::writeattribute(group_id , "x_pts"      , x_pts);
	ionsim::writeattribute(group_id , "y_pts"      , y_pts);
	ionsim::writeattribute(group_id , "x_edge_mag" , x_edge_mag);
	ionsim::writeattribute(group_id , "y_edge_mag" , y_edge_mag);
	
	// ==================================
	// Close file
	// ==================================
	/* hdf5_cleanup(x_dataspace_id, y_dataspace_id, x_dataset_id, y_dataset_id, file_id); */
	H5Gclose(step_group_id);
	H5Gclose(group_id);
	H5Sclose(x_dataspace_id);
	H5Sclose(y_dataspace_id);
	H5Dclose(x_dataset_id);
	H5Dclose(y_dataset_id);
	H5Sclose(x_memspace_id);
	H5Sclose(y_memspace_id);
	H5Fclose(file_id);

	
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

int Field::recv_field_others()
{
	double xbuf[n_pts];
	double ybuf[n_pts];
	int p = MPI::COMM_WORLD.Get_size();
	int my_id = MPI::COMM_WORLD.Get_rank();

	for (int id=0; id < p; id++)
	{
		if (id != my_id) {
			MPI::COMM_WORLD.Recv(xbuf, n_pts, MPI::DOUBLE, id, ionsim::TAG_FIELD);
			MPI::COMM_WORLD.Recv(ybuf, n_pts, MPI::DOUBLE, id, ionsim::TAG_FIELD);

			for (int j=0; j < n_pts; j++)
			{
				x_data[j] += xbuf[j];
				y_data[j] += ybuf[j];
			}
		}
	}

	splines_valid = false;

	return 0;
}

int Field::recv_field(int sender_id)
{
	double xbuf[n_pts];
	double ybuf[n_pts];
	int p = MPI::COMM_WORLD.Get_size();
	int my_id = MPI::COMM_WORLD.Get_rank();

	MPI::COMM_WORLD.Recv(xbuf, n_pts, MPI::DOUBLE, sender_id, ionsim::TAG_FIELD);
	MPI::COMM_WORLD.Recv(ybuf, n_pts, MPI::DOUBLE, sender_id, ionsim::TAG_FIELD);

	for (int j=0; j < n_pts; j++)
	{
		x_data[j] += xbuf[j];
		y_data[j] += ybuf[j];
	}

	splines_valid = false;

	return 0;
}

int Field::send_field(int dest_id)
{
	MPI::COMM_WORLD.Send(x_data, n_pts, MPI::DOUBLE, dest_id, ionsim::TAG_FIELD);
	MPI::COMM_WORLD.Send(y_data, n_pts, MPI::DOUBLE, dest_id, ionsim::TAG_FIELD);

	return 0;
}

int Field::get_interp(Field &field)
{
	for (int i=0; i < field.x_pts; i++)
	{
		for (int j=0; j < field.y_pts; j++)
		{
			field.Ex_ind(i, j) = Ex(field.i_to_x(i), field.j_to_y(j));
			field.Ey_ind(i, j) = Ey(field.i_to_x(i), field.j_to_y(j));
		}
	}

	return 0;
}
