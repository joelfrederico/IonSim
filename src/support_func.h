#ifndef __SUPPORT_FUNC_H_INCLUDED__
#define __SUPPORT_FUNC_H_INCLUDED__

#include <mpi.h>
#include <string>
/* #include "parts.h" */
#include <hdf5.h>
#include "field_data.h"
#include "field_comm.h"
#include <numeric>

namespace ionsim
{
	// ==================================
	// Methods
	// ==================================
	double GeV2gamma(double GeV);
	double gamma2GeV(double gamma);
	double gaussian(double mean, double sigma);

	double sr(double emit_n, double E, double n_p_cgs, double m_ion_amu);
	double nb_0(double q_tot, double sz, double sr);
	double nb_0(double q_tot, double sz, double emit_n, double E, double n_p_cgs, double m_ion_amu);

	template<typename T>
	T mean(std::vector<T> vec)
	{
		return std::accumulate(vec.begin(), vec.end(), 0) / vec.size();
	}

	template<typename T>
	T std(const std::vector<T> vec, const T &mean)
	{
		std::vector<T> temp(vec.size());
		std::transform(vec.begin(), vec.end(), temp.begin(), [&mean](const T &val) { return val - mean;});
		return std::inner_product(vec.begin(), vec.end(), vec.begin(), 0) / vec.size();
	}

	/* template<typename T> */
	/* T col_major(const T i, const T j, const ptrdiff_t N0) */
	/* { */
	/* 	return i + N0*j; */
	/* } */

	template<typename Tclass>
	auto row_major(const std::vector<Tclass> _x_pts, const typename decltype(_x_pts)::value_type i) -> typename decltype(_x_pts)::value_type
	{
		return i;
	}

	template<typename Tclass>
	auto row_major(const std::vector<Tclass> _x_pts) -> typename decltype(_x_pts)::value_type
	{
		return 1;
	}

	template<typename Tclass, typename... T2>
	auto row_major(const std::vector<Tclass> _x_pts, const typename decltype(_x_pts)::value_type i, const T2 ... rest) -> typename decltype(_x_pts)::value_type
	{
		auto size = sizeof...(rest);
		typename decltype(_x_pts)::size_type i_dim;
		typename decltype(_x_pts)::value_type f;
	
		i_dim = _x_pts.size() - size - 1;
		/* ind = size; */
	
		if (i >= _x_pts[i_dim]) throw std::runtime_error("Index requested is too large in a dimension.");
	
		// Column-major
		/* return i + _x_pts[ind] * _index(rest...); */

		// Row-major
		f = 1;
		for (decltype(size) ii=i_dim+1; ii < _x_pts.size(); ii++)
		{
			f *= _x_pts[ii];
		}
		return i * f + row_major(_x_pts, rest...);
	}

	
	/* template<typename T> */
	/* T row_major(const T i, const T j, const ptrdiff_t N1) */
	/* { */
	/* 	return i*N1 + j; */
	/* } */

	template<typename T>
	MPI_Datatype convert_typeid_to_mpi()
	{
		MPI_Datatype mpi_type;
		if (typeid(T) == typeid(long double))
		{
			mpi_type = MPI_LONG_DOUBLE;
		} else if (typeid(T) == typeid(unsigned int)) {
			mpi_type = MPI_UNSIGNED;
		} else if (typeid(T) == typeid(int)) {
			mpi_type = MPI_INT;
		} else {
			throw std::runtime_error("Type not handled!");
		}
		return mpi_type;
	}
}

#endif
