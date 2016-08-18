#ifndef __SCALAR_DATA_H_INCLUDED__
#define __SCALAR_DATA_H_INCLUDED__

#include "consts.h"
#include "simparams.h"

template<typename Tclass>
class ScalarData
{
	private:
		// ==================================
		// Private data members
		// ==================================
		// Initialization parameters
		std::vector<std::vector<long double>::size_type> _x_pts;
		std::vector<long double> 	                 _edge_mag;

		// Metadata parameters
		std::vector<double>      	      mid;
		std::vector<long double> 	      _dxdi;
		std::vector<std::vector<long double>> _grid;

		// Data storage
		std::vector<Tclass>      	      data;

		// ==================================
		// Private methods
		// ==================================
		bool _samedim(const ScalarData<Tclass> &rhs) const;

		int _init(
				const std::vector<typename std::vector<Tclass>::size_type> x_pts,
				const std::vector<long double> edge_mag
			);


	public:
		// ==================================
		// Constructors, Destructor
		// ==================================
		ScalarData();
		ScalarData(const SimParams &simparams);
		ScalarData(const decltype(_x_pts) x_pts, const decltype(_edge_mag) edge_mag);

		// ==================================
		// Calculated data
		// ==================================
		auto dxdi     	  (const decltype(_dxdi)::size_type i)     const -> typename decltype(_dxdi)::value_type;
		auto x_pts    	  (const decltype(_x_pts)::size_type i)    const -> typename decltype(_x_pts)::value_type;
		auto x_pts_vec	  ()                                       const -> decltype(_x_pts);
		auto n_pts    	  () 				           const -> typename decltype(data)::size_type;
		auto edge_mag     (const decltype(_edge_mag)::size_type i) const -> typename decltype(_edge_mag)::value_type;
		auto edge_mag_vec ()                                       const -> decltype(_edge_mag);
		auto grid     	  (const decltype(_grid)::size_type i)     const -> typename decltype(_grid)::value_type;
		auto vdata    	  ()  				       const -> decltype(data);

		// ==================================
		// Methods
		// ==================================
		void resize(const decltype(_x_pts) x_pts, const decltype(_edge_mag) edge_mag);

		template<typename T, typename... T2>
		auto _index(const T i, const T2 ... rest) const -> typename decltype(data)::size_type
		{
			auto size = sizeof...(rest);
			typename decltype(_x_pts)::size_type ind;
		
			ind = _x_pts.size() - size - 1;
		
			if (i >= _x_pts[ind]) throw std::runtime_error("Index requested is too large in a dimension.");
		
			return i + _x_pts[ind] * _index(rest...);
		}

		auto _index(const typename decltype(data)::size_type i) const -> typename decltype(data)::size_type
		{
			return i;
		}

		template<typename... T2>
		Tclass &ind(const T2 ...rest)
		{
			auto ind = _index(rest...);
			return data[ind];
		}

		Tclass &ind(typename decltype(data)::size_type i);

		int lt_x_ind_e(const typename decltype(_x_pts)::size_type i, const long double x, decltype(_x_pts)::value_type &ind) const;

		// ==================================
		// Operators
		// ==================================
		/* ScalarData<Tclass> &operator=(const ScalarData<Tclass> &rhs); */
		ScalarData<Tclass> &operator+=(const ScalarData<Tclass> &rhs);
		ScalarData<Tclass> &operator-=(const ScalarData<Tclass> &rhs);
		ScalarData<Tclass> &operator=(const Tclass rhs);

		ScalarData<Tclass> &operator*=(const Tclass rhs);
		/* { */
		/* 	for (int i=0; i < n_pts(); i++) */
		/* 	{ */
		/* 		this->data[i] *= rhs; */
		/* 	} */
		/* 	return *this; */
		/* } */

		const ScalarData<Tclass> operator+(const ScalarData<Tclass> &rhs);
		const ScalarData<Tclass> operator-(const ScalarData<Tclass> &rhs);
		ScalarData<Tclass> operator-() const;

		template <class T>
		const ScalarData<Tclass> operator*(const T rhs)
		{
			return ScalarData<Tclass>(*this) *= rhs;
		}
};

#endif
