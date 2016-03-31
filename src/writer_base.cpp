#include "writer_base.h"
#include "field_data.h"
#include "support_func.h"

WriterBase::WriterBase(const std::string &filename)
{
}

WriterBase::~WriterBase()
{
	H5Fclose(file_id);
}


hid_t WriterBase::group_access(hid_t &loc_id, std::string const &group)
{
	return group_access(loc_id, group.c_str());
}

hid_t WriterBase::group_access(hid_t &loc_id, const char *group)
{
	herr_t status;
	hid_t group_id;
	H5G_info_t objinfo;

	// ==================================
	// Access or create a new group
	// ==================================
	H5Eset_auto(NULL, NULL, NULL);
	group_id = H5Gopen(loc_id, group, H5P_DEFAULT);
	if (group_id < 0)
	{
		group_id = H5Gcreate(loc_id, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gget_info_by_name(loc_id, group, &objinfo, H5P_DEFAULT);
	}
	return group_id;
}

hid_t WriterBase::group_step_access(hid_t &loc_id, long step)
{
	char buf[10];
	sprintf(buf, "Step_%04ld", step);
	return group_access(loc_id, buf);
}

hid_t WriterBase::dataset_create(hid_t &group_id, hid_t &dataspace_id, hsize_t count[2], std::string const &dataset)
{
	hid_t dataset_id;

	dataspace_id = H5Screate_simple(2, count, NULL);

	// ==================================
	// Create dataset
	// ==================================
	dataset_id = H5Dcreate(group_id, dataset.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	return dataset_id;
}

int WriterBase::dump_serial(const Field_Data &field, long step)
{
	/* ionsim::dump_serial(filename, step, group, *this); */
	// ==================================
	// Initialize all variables
	// ==================================
	double *xbuf;
	double *ybuf;
	double *zbuf;
	std::string x_dataset, y_dataset, z_dataset;
	long x_len = field.x_pts;
	long y_len = field.y_pts;
	long z_len = field.z_pts;

	int n_write = ionsim::MAX_N_WRITE;

	hid_t step_group_id, group_id;
	hid_t x_dataspace_id, y_dataspace_id;
	hid_t x_dataset_id, y_dataset_id;
	hid_t x_memspace_id, y_memspace_id;
	hsize_t count[3];
	hsize_t offset[3];
	std::string temp_str;
	herr_t status;
	
	// ==================================
	// Access or create a new group
	// ==================================
	step_group_id    = group_step_access(file_id, step);

	group_id         = group_access(step_group_id, "field");

	// ==================================
	// Create dataset
	// ==================================
	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;

	x_dataset = "Ex";
	x_dataset_id = dataset_create(group_id, x_dataspace_id, count, x_dataset); // Updates dataspace_id
	x_memspace_id = H5Screate_simple(2, count, NULL);

	count[0] = x_len;
	count[1] = y_len;
	y_dataset = "Ey";
	y_dataset_id = dataset_create(group_id, y_dataspace_id, count, y_dataset); // Updates dataspace_id
	y_memspace_id = H5Screate_simple(2, count, NULL);

	// ==================================
	// Write dataset rows
	// ==================================
	xbuf = new double[x_len*y_len*z_len];
	ybuf = new double[x_len*y_len*z_len];
	zbuf = new double[x_len*y_len*z_len];
	for (int i=0; i < x_len; i++)
	{
		for (int j=0; j < y_len; j++)
		{
			for (int k=0; k < z_len; k++)
			{
				xbuf[i + x_len*(j + k*y_len)] = field.Ex_ind(i, j, k);
				ybuf[i + x_len*(j + k*y_len)] = field.Ey_ind(i, j, k);
				zbuf[i + x_len*(j + k*y_len)] = field.Ez_ind(i, j, k);
			}
		}
	}
	H5Dwrite(x_dataset_id, H5T_NATIVE_DOUBLE, x_memspace_id, x_dataspace_id, H5P_DEFAULT, xbuf);
	H5Dwrite(y_dataset_id, H5T_NATIVE_DOUBLE, y_memspace_id, y_dataspace_id, H5P_DEFAULT, ybuf);

	delete [] xbuf;
	delete [] ybuf;

	writeattribute(group_id , "n_pts"      , field.n_pts);
	/* writeattribute(group_id , "dxdi"       , field.dxdi); */
	/* writeattribute(group_id , "dydj"       , field.dydj); */
	/* writeattribute(group_id , "mid_i"      , field.mid_i); */
	/* writeattribute(group_id , "mid_j"      , field.mid_j); */
	writeattribute(group_id , "x_pts"      , field.x_pts);
	writeattribute(group_id , "y_pts"      , field.y_pts);
	writeattribute(group_id , "z_pts"      , field.z_pts);
	writeattribute(group_id , "x_edge_mag" , field.x_edge_mag);
	writeattribute(group_id , "y_edge_mag" , field.y_edge_mag);
	writeattribute(group_id , "z_edge_mag" , field.z_edge_mag);
	
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
	
	return 0;
}

int WriterBase::write_attributes_parallel(const SimParams &simparams) const
{

	writeattribute(file_id, "n_e"       , simparams.n_e      );
	writeattribute(file_id, "n_ions"    , simparams.n_ions   );
	writeattribute(file_id, "q_tot"     , simparams.q_tot    );
	writeattribute(file_id, "radius"    , simparams.radius   );
	writeattribute(file_id, "length"    , simparams.length   );
	writeattribute(file_id, "E"         , simparams.E        );
	writeattribute(file_id, "emit_n"    , simparams.emit_n   );
	writeattribute(file_id, "n_p_cgs"   , simparams.n_p_cgs  );
	writeattribute(file_id, "m_ion_amu" , simparams.m_ion_amu);
	writeattribute(file_id, "sz"        , simparams.sz       );
	writeattribute(file_id, "sdelta"    , simparams.sdelta   );
	writeattribute(file_id, "dt"        , simparams.dt       );
	writeattribute(file_id, "n_steps"   , simparams.n_steps  );
	writeattribute(file_id, "n_field_x" , simparams.n_field_x);
	writeattribute(file_id, "n_field_y" , simparams.n_field_y);
	writeattribute(file_id, "n_field_z" , simparams.n_field_z);
	H5Fclose(file_id);

	return 0;
}

int WriterBase::dump_serial(std::string const &filename, long step, std::string const &group, const Field_Data &field)
{
	// ==================================
	// Initialize all variables
	// ==================================
	double *xbuf;
	double *ybuf;
	double *zbuf;
	std::string x_dataset, y_dataset, z_dataset;
	long x_len = field.x_pts;
	long y_len = field.y_pts;
	long z_len = field.z_pts;

	int n_write = ionsim::MAX_N_WRITE;

	hid_t step_group_id, group_id;
	hid_t x_dataspace_id, y_dataspace_id, z_dataspace_id;
	hid_t x_dataset_id, y_dataset_id, z_dataset_id;
	hid_t x_memspace_id, y_memspace_id, z_memspace_id;
	hsize_t count[3];
	hsize_t offset[3];
	std::string temp_str;
	herr_t status;

	// ==================================
	// Access or create a new group
	// ==================================
	step_group_id    = group_step_access(file_id, step);

	group_id         = group_access(step_group_id, group);

	// ==================================
	// Create dataset
	// ==================================
	// The total array of particles is (n_pts*num_processes, 6) in size
	// Remember that each slave has n_pts of different particles!
	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;
	x_dataset = "Ex";
	x_dataset_id = dataset_create(group_id, x_dataspace_id, count, x_dataset); // Updates dataspace_id
	x_memspace_id = H5Screate_simple(2, count, NULL);

	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;
	y_dataset = "Ey";
	y_dataset_id = dataset_create(group_id, y_dataspace_id, count, y_dataset); // Updates dataspace_id
	y_memspace_id = H5Screate_simple(2, count, NULL);

	count[0] = x_len;
	count[1] = y_len;
	count[2] = z_len;
	z_dataset = "Ez";
	z_dataset_id = dataset_create(group_id, z_dataspace_id, count, z_dataset); // Updates dataspace_id
	z_memspace_id = H5Screate_simple(2, count, NULL);

	// ==================================
	// Write dataset rows
	// ==================================
	xbuf = new double[x_len*y_len*z_len];
	ybuf = new double[x_len*y_len*z_len];
	zbuf = new double[x_len*y_len*z_len];

	for (int i=0; i < x_len; i++)
	{
		for (int j=0; j < y_len; j++)
		{
			for (int k=0; k < z_len; k++)
			{
				xbuf[i + x_len*(j + k*y_len)] = field.Ex_ind(i, j, k);
				ybuf[i + x_len*(j + k*y_len)] = field.Ey_ind(i, j, k);
				zbuf[i + x_len*(j + k*y_len)] = field.Ez_ind(i, j, k);
			}
		}
	}
	H5Dwrite(x_dataset_id, H5T_NATIVE_DOUBLE, x_memspace_id, x_dataspace_id, H5P_DEFAULT, xbuf);
	H5Dwrite(y_dataset_id, H5T_NATIVE_DOUBLE, y_memspace_id, y_dataspace_id, H5P_DEFAULT, ybuf);
	H5Dwrite(z_dataset_id, H5T_NATIVE_DOUBLE, z_memspace_id, z_dataspace_id, H5P_DEFAULT, zbuf);

	delete [] xbuf;
	delete [] ybuf;

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

	return 0;
}

