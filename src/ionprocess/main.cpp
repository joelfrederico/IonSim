/* #include "config.h" */
#include "main.h"
#include <gflags/gflags.h> // #include <google/gflags.h>
#include <ionprocess.h>
#include "../hdf5_classes.h"

DEFINE_bool(verbose, true, "Verbose");
DEFINE_int32(step, 0, "Electron Step");
DEFINE_int32(xbins, 101, "Number of bins in x");
/* DEFINE_string(file, "config.h5", "Data file"); */



int main(int argc, char **argv)
{
	std::string filename, basename, procname;
	hsize_t dims[2];
	hsize_t memdims[1];

	// ==============================
	// Parse flags
	// ==============================
	gflags::ParseCommandLineFlags(&argc, &argv, true);
	int xbins = FLAGS_xbins;

	// ==============================
	// Get filename
	// ==============================
	if (argc < 2)
	{
		std::cout << "Not enough arguments!" << std::endl;
		return 0;
	} else {
		filename = argv[1];
	}

	// ==============================
	// Get extract base, proc names
	// ==============================
	basename = filename.substr(0, filename.find_last_of("."));
	procname = basename + "_processed.h5";

	unsigned long long *hist;
	long histsize;
	int n_field_z;

	makehist(filename, xbins, 0, hist, histsize, n_field_z);
	/* std::cout << "Hist[0]: " << hist[0] << std::endl; */
	/* std::cout << "Hist: " << hist << std::endl; */
	/* std::cout << "histsize: " << histsize << std::endl; */
	/* std::cout << "n_field_z: " << n_field_z << std::endl; */

	// ==============================
	// Write to file
	// ==============================
	DataspaceCreate *dataspace, *memspace;
	FileCreate output(procname);
	dims[0] = xbins;
	dims[1] = n_field_z;
	memdims[0] = xbins*n_field_z;
	DatasetAccess hist_data(output.file_id, "hist", 2, dims);
	memspace = new DataspaceCreate(1, memdims);

	H5Dwrite(hist_data.dataset_id, H5T_NATIVE_ULLONG, memspace->dataspace_id, hist_data.dataspace_id, H5P_DEFAULT, hist);

	delete memspace;
	delete[] hist;
	return 0;
}
