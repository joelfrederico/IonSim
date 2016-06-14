#include "../hdf5_classes.h"
#include <cmath>
#include "ionprocess.h"
#include <vector>
#include <string>

std::string _getion(unsigned int step)
{
	std::string out;
	char buf[10];
	sprintf(buf, "ions_%04u", step);
	out = buf;
	return out;
}

int index(int i, int j, int ncols)
{
	return i + j * ncols;
}

int index(int i, int j, int k, int ncols, int nrows)
{
	return i + j * ncols + k*ncols*nrows;
}

int makehist(std::string filename, int xbins, int step, std::vector<unsigned long long> &hist, long &histsize, int &n_field_z)
{
	// ==============================
	// Initialize Vars 
	// ==============================
	DatasetOpen *ions;
	DataspaceCreate *dataspace, *memspace;
	double *pts;
	bool *pts_bool;

	double dx = 1e-5;
	double delx;

	herr_t status;
	hsize_t dims[2];
	hsize_t memdims[1];
	int ind;
	int n_filtered;
	double z_end;
	double x, xp, y, yp, z, zp;
	double x_abs_max;
	double *count, *z_array;

	// ==============================
	// Open files for reading/writing
	// ==============================
	FileOpen data(filename);

	// ==============================
	// Find ions and number of them
	// ==============================
	GroupStepAccess stepgroup(data.file_id, step);
	GroupAccess ionstep(stepgroup.group_id, "ions_steps");

	AttributeOpen n_field_z_attr(data.file_id, "n_field_z");
	n_field_z = n_field_z_attr.read();
	AttributeOpen z_end_attr(data.file_id, "z_end");
	z_end = z_end_attr.read();

	// ==============================
	// Loop over ions to find range
	// ==============================
	ions      = new DatasetOpen(ionstep.group_id, _getion(0));
	H5Sget_simple_extent_dims(ions->dataspace_id, dims, NULL);
	memdims[0] = dims[0]*dims[1];
	delete ions;

	x_abs_max = 0;
	n_filtered = 0;
	for (int i=0; i < n_field_z ; i++)
	{
		ions      = new DatasetOpen(ionstep.group_id, _getion(i));
		H5Sget_simple_extent_dims(ions->dataspace_id, dims, NULL);

		pts = new double[memdims[0]];

		dataspace = new DataspaceCreate(2, dims);
		memspace  = new DataspaceCreate(1, memdims);

		status = H5Dread(ions->dataset_id, H5T_NATIVE_DOUBLE, memspace->dataspace_id, dataspace->dataspace_id, H5P_DEFAULT, &(*pts));

		for (int j=0; j < dims[0]; j++)
		{
			x = pts[index(0, j, 6)];
			if ((std::abs(x) <= dx) && (std::abs(x) > x_abs_max))
			{
				n_filtered++;
				/* std::cout << i << std::endl; */
				x_abs_max = std::abs(x);
			}
		}

		delete [] pts;
		delete memspace;
		delete dataspace;
		delete ions;
	}
	std::cout << "Max: " << x_abs_max << std::endl;

	// ==============================
	// Histogram
	// ==============================
	histsize = xbins*n_field_z;
	hist.reserve(histsize);
	delx = 2*x_abs_max/xbins;

	for (int i=0; i < n_field_z ; i++)
	{
		ions      = new DatasetOpen(ionstep.group_id, _getion(i));
		H5Sget_simple_extent_dims(ions->dataspace_id, dims, NULL);

		pts = new double[memdims[0]];

		dataspace = new DataspaceCreate(2, dims);
		memspace  = new DataspaceCreate(1, memdims);

		status = H5Dread(ions->dataset_id, H5T_NATIVE_DOUBLE, memspace->dataspace_id, dataspace->dataspace_id, H5P_DEFAULT, &(*pts));

		for (int j=0; j < dims[0]; j++)
		{
			x = pts[index(0, j, 6)];
			if (std::abs(x) <= dx)
			{
				ind = floor((x+x_abs_max) / delx);
				hist[index(i, ind, n_field_z)]++;
			}
		}

		delete [] pts;
		delete memspace;
		delete dataspace;
		delete ions;
	}

	std::cout << "xbins: " << xbins << std::endl;
	std::cout << "n_field_z: " << n_field_z << std::endl;
	std::cout << "step: " << step << std::endl;
	std::cout << "Hist[0]: " << hist[0] << std::endl;

	unsigned long long test = 0;
	std::cout << "Test: " << test - 1 << std::endl;

	return 0;
}
