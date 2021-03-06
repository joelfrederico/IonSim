#include "config.h"
#include "writer_base.h"
#include "field_data.h"
#include "support_func.h"
#include <ionsim.h>

WriterBase::WriterBase(const std::string &filename)
{
	_filename = filename;
	writer_type = WRITER_NULL;
}

WriterBase::~WriterBase()
{
	herr_t status;
	ssize_t num_open;
	hid_t *obj_id_list;

	status = H5Fclose(file_id);
	if (status == -1)
	{
		num_open = H5Fget_obj_count(file_id, H5F_OBJ_ALL);
		obj_id_list = new hid_t[num_open];
		H5Fget_obj_ids(file_id, H5F_OBJ_ALL, num_open, obj_id_list);

		std::cout << "ERROR: (" << status << "). COULD NOT CLOSE FILE." << std::endl;
		std::cout << "Open obj count: " << num_open << std::endl;

		for (int i=0; i < num_open; i++)
		{
			std::cout << "Obj id: " << obj_id_list[i] << std::endl;
		}

		delete obj_id_list;
	}

}

hid_t WriterBase::dataset_create(hid_t &group_id, hid_t &dataspace_id, int rank, hsize_t *count, const std::string &dataset)
{
	hid_t dataset_id;

	dataspace_id = H5Screate_simple(rank, count, NULL);

	// ==================================
	// Create dataset
	// ==================================
	dataset_id = H5Dcreate(group_id, dataset.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	return dataset_id;
}

int WriterBase::write_attributes(const SimParams &simparams) const
{
	AttributeCreate version          ( file_id , "version"          , VERSION                    );
	AttributeCreate n_e              ( file_id , "n_e"              , simparams.n_e              );
	AttributeCreate n_e_node         ( file_id , "n_e_node"         , simparams.n_e_node()       );
	AttributeCreate n_ions           ( file_id , "n_ions"           , simparams.n_ions           );
	AttributeCreate q_tot            ( file_id , "q_tot"            , simparams.q_tot            );
	AttributeCreate qpp_e            ( file_id , "qpp_e"            , simparams.qpp_e()          );
	AttributeCreate radius           ( file_id , "radius"           , simparams.radius           );
	AttributeCreate E                ( file_id , "E"                , simparams.E                );
	AttributeCreate emit_n           ( file_id , "emit_n"           , simparams.emit_n           );
	AttributeCreate n_p_cgs          ( file_id , "n_p_cgs"          , simparams.n_p_cgs          );
	AttributeCreate m_ion_amu        ( file_id , "m_ion_amu"        , simparams.m_ion_amu        );
	AttributeCreate sz               ( file_id , "sz"               , simparams.sz               );
	AttributeCreate sdelta           ( file_id , "sdelta"           , simparams.sdelta           );
	AttributeCreate dt               ( file_id , "dt"               , simparams.dt()             );
	AttributeCreate n_steps          ( file_id , "n_steps"          , simparams.n_steps          );
	AttributeCreate n_field_x        ( file_id , "n_field_x"        , simparams.n_field_x        );
	AttributeCreate n_field_y        ( file_id , "n_field_y"        , simparams.n_field_y        );
	AttributeCreate n_field_z        ( file_id , "n_field_z"        , simparams.n_field_z        );
	AttributeCreate field_trans_wind ( file_id , "field_trans_wind" , simparams.field_trans_wind );
	AttributeCreate z_end            ( file_id , "z_end"            , simparams.z_end            );

	AttributeCreate pushmethod_int   ( file_id , "pushmethod_int"   , simparams.pushmethod       );
	PushMethod pushmeth(simparams.pushmethod);
	AttributeCreate pushmethod_name  ( file_id , "pushmethod_name"  , pushmeth.name.c_str()      );

	AttributeCreate zdist_int        ( file_id , "zdist_int"        , simparams.zdist            );
	zDist zdist(simparams.zdist);
	AttributeCreate zdist_name  ( file_id , "zdist_name"  , zdist.name.c_str()      );


	return 0;
}

