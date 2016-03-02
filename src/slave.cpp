#include "mpi.h"
#include "classes.h"

int slave(int &p, int &id, MPI::Intracomm &slave_comm_id)
{

	// ==============================
	// Receive starting gun
	// ==============================
	int whoami;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(&whoami, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
	printf("Process %d says: **DUM DUM DUM DUM**", id);

	// ==============================
	// Generate beam
	// ==============================
	int root = 0;
	int n_e;
	int n_ion;
	double q_tot;
	double x_window;
	double y_window;
	double E;
	double emit_n;
	double n_p_cgs;
	double m_ion_amu;
	double sz;
	double sdelta;

	MPI::COMM_WORLD.Bcast(&n_e, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&n_ion, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&x_window, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&y_window, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&E, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&emit_n, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&n_p_cgs, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&m_ion_amu, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&sz, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&sdelta, 1, MPI::DOUBLE, root);

	double z_cov[2][2] = {{sz, 0}, {0, sdelta}};

	Emit emit(emit_n, E, true);
	Plasma plas(n_p_cgs, m_ion_amu);
	Match mat(plas, E, emit);
	Beam x_beam(mat.beta(), mat.alpha(), emit);
	Beam y_beam(mat.beta(), mat.alpha(), emit);
	Ions ions(plas, n_ion, x_window, y_window);

	printf("Id: %d, n_e: %d, n_ion: %d, x_win: %.3e, y_win: %.3e, E: %.3e\n", id, n_e, n_ion, x_window, y_window, E);

	Ebeam ebeam(n_e, q_tot, E, x_beam, y_beam, z_cov);
	ebeam.dump("output.h5", slave_comm_id);

	return 0;
}
