#include "config.h"
#include "main.h"
#include "../hdf5_classes.h"
#include <gflags/gflags.h> // #include <google/gflags.h>

DEFINE_bool(verbose, true, "Verbose");
DEFINE_int32(step, 0, "Electron Step");
/* DEFINE_string(file, "config.h5", "Data file"); */


std::string _getion(unsigned int step)
{
	std::string out;
	char buf[10];
	sprintf(buf, "ions_%04u", step);
	out = buf;
	return out;
}

int main(int argc, char **argv)
{
	// ==============================
	// Initialize Vars 
	// ==============================
	std::string filename;
	std::string basename;
	std::string procname;

	// ==============================
	// Parse flags
	// ==============================
	gflags::ParseCommandLineFlags(&argc, &argv, true);

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

	// ==============================
	// Open files for reading/writing
	// ==============================
	FileOpen data(filename);
	FileCreate output(procname);

	// ==============================
	// Find ions and number of them
	// ==============================
	GroupStepAccess step(data.file_id, FLAGS_step);
	GroupAccess ionstep(step.group_id, "ions_steps");
	AttributeOpen n_ions(data.file_id, "n_field_z");

	std::cout << "N_ions: " << n_ions.read() << std::endl;
	
	// ==============================
	// Loop over ions
	// ==============================
	DatasetOpen *ions;
	ions = new DatasetOpen(ionstep.group_id, _getion(0));
	delete ions;

}
